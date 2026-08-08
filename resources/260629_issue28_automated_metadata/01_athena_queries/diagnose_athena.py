import os
import json
import subprocess
import time

def run_command(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\nError: {result.stderr}")
    return result.stdout.strip()

def check_query_status(execution_id):
    cmd = f"aws athena get-query-execution --query-execution-id {execution_id}"
    output = run_command(cmd)
    data = json.loads(output)
    return data["QueryExecution"]["Status"]["State"], data["QueryExecution"]["Status"].get("StateChangeReason", "")

def run_athena_query(query_string):
    output_location = "s3://petabite/automated-metadata/athena-results/"
    payload = {
        "QueryString": query_string,
        "QueryExecutionContext": {"Database": "default"},
        "ResultConfiguration": {"OutputLocation": output_location}
    }
    temp_file = "temp_diag_query.json"
    with open(temp_file, "w") as f:
        json.dump(payload, f)
    try:
        cmd = f"aws athena start-query-execution --cli-input-json file://{temp_file}"
        output = run_command(cmd)
        execution_id = json.loads(output)["QueryExecutionId"]
        while True:
            state, reason = check_query_status(execution_id)
            if state == "SUCCEEDED":
                # Get query results
                get_results_cmd = f"aws athena get-query-results --query-execution-id {execution_id}"
                results_output = run_command(get_results_cmd)
                return json.loads(results_output)
            elif state in ["FAILED", "CANCELLED"]:
                raise RuntimeError(f"Query {execution_id} {state}. Reason: {reason}")
            time.sleep(2)
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)

def extract_single_value(results):
    try:
        return results["ResultSet"]["Rows"][1]["Data"][0]["VarCharValue"]
    except (IndexError, KeyError):
        return "0"

def main():
    print("Starting Athena Diagnostics for prefix SAMN14...")
    
    # Check 1: Count in default.petadex_biosamples for SAMN14
    q1 = "SELECT COUNT(*) FROM default.petadex_biosamples WHERE biosample LIKE 'SAMN14%';"
    try:
        r1 = run_athena_query(q1)
        val1 = extract_single_value(r1)
        print(f"[Check 1] Rows in default.petadex_biosamples starting with SAMN14: {val1}")
    except Exception as e:
        print(f"[Check 1] Failed: {e}")

    # Check 2: Count in sra.metadata for SAMN14
    q2 = "SELECT COUNT(*) FROM sra.metadata WHERE biosample LIKE 'SAMN14%';"
    try:
        r2 = run_athena_query(q2)
        val2 = extract_single_value(r2)
        print(f"[Check 2] Rows in sra.metadata starting with SAMN14: {val2}")
    except Exception as e:
        print(f"[Check 2] Failed: {e}")

    # Check 3: Joined count (without tax_analysis) for SAMN14
    q3 = """
    SELECT COUNT(*) 
    FROM sra.metadata meta
    JOIN default.petadex_biosamples pb ON meta.biosample = pb.biosample
    WHERE meta.biosample LIKE 'SAMN14%';
    """
    try:
        r3 = run_athena_query(q3)
        val3 = extract_single_value(r3)
        print(f"[Check 3] Joined rows (meta + pb) for SAMN14: {val3}")
    except Exception as e:
        print(f"[Check 3] Failed: {e}")

    # Check 4: Check if SRA Tax Analysis contains matches
    q4_sample = "SELECT acc FROM sra.metadata WHERE biosample LIKE 'SAMN14%' LIMIT 5;"
    try:
        r4_sample = run_athena_query(q4_sample)
        rows = r4_sample["ResultSet"]["Rows"]
        if len(rows) > 1:
            run_ids = [row["Data"][0]["VarCharValue"] for row in rows[1:]]
            print(f"[Check 4] Sample Run IDs found in sra.metadata: {run_ids}")
            run_ids_str = ", ".join(f"'{rid}'" for rid in run_ids)
            
            q4 = f"SELECT COUNT(*) FROM sra_tax_analysis_tool.tax_analysis WHERE acc IN ({run_ids_str});"
            r4 = run_athena_query(q4)
            val4 = extract_single_value(r4)
            print(f"[Check 4] Tax analysis rows for these runs: {val4}")
        else:
            print("[Check 4] No runs found in sra.metadata for SAMN14 to test tax_analysis.")
    except Exception as e:
        print(f"[Check 4] Failed: {e}")

if __name__ == '__main__':
    main()
