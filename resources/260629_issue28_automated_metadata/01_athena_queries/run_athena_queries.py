import os
import json
import subprocess
import time

def run_command(cmd):
    """Run a shell command and return stdout as string, or raise error on failure."""
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {cmd}\nError: {result.stderr}")
    return result.stdout.strip()

def check_query_status(execution_id):
    """Check the status of an Athena query execution."""
    cmd = f"aws athena get-query-execution --query-execution-id {execution_id}"
    output = run_command(cmd)
    data = json.loads(output)
    return data["QueryExecution"]["Status"]["State"], data["QueryExecution"]["Status"].get("StateChangeReason", "")

def run_athena_query(query_string, output_location):
    """Submit a query to Athena, wait for completion, and return the execution ID."""
    payload = {
        "QueryString": query_string,
        "QueryExecutionContext": {"Database": "default"},
        "ResultConfiguration": {"OutputLocation": output_location}
    }
    
    temp_file = "temp_query_payload.json"
    with open(temp_file, "w") as f:
        json.dump(payload, f)
        
    try:
        cmd = f"aws athena start-query-execution --cli-input-json file://{temp_file}"
        output = run_command(cmd)
        execution_id = json.loads(output)["QueryExecutionId"]
        
        while True:
            state, reason = check_query_status(execution_id)
            if state == "SUCCEEDED":
                return execution_id
            elif state in ["FAILED", "CANCELLED"]:
                raise RuntimeError(f"Query {execution_id} {state}. Reason: {reason}")
            time.sleep(2)
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)

def main():
    queries_dir = "resources/260629_issue28_automated_metadata/queries"
    output_location = "s3://petabite/automated-metadata/athena-results/"
    
    if not os.path.exists(queries_dir):
        print(f"Queries directory {queries_dir} not found. Please run query generator first.")
        return
        
    sql_files = [f for f in os.listdir(queries_dir) if f.endswith(".sql")]
    print(f"Found {len(sql_files)} SQL queries in {queries_dir}.")
    
    completed_log = "resources/260629_issue28_automated_metadata/completed_queries.txt"
    completed = set()
    if os.path.exists(completed_log):
        with open(completed_log, "r") as f:
            completed = set(line.strip() for line in f if line.strip())
            
    print(f"Skipping {len(completed)} queries that are already completed.")
    
    for i, file_name in enumerate(sorted(sql_files)):
        prefix = file_name.replace(".sql", "")
        if prefix in completed:
            continue
            
        file_path = os.path.join(queries_dir, file_name)
        with open(file_path, "r") as f:
            query_string = f.read()
            
        print(f"[{i+1}/{len(sql_files)}] Running query for prefix {prefix}...")
        start_time = time.time()
        try:
            exec_id = run_athena_query(query_string, output_location)
            elapsed = time.time() - start_time
            print(f"-> Succeeded in {elapsed:.2f}s. ID: {exec_id}")
            
            with open(completed_log, "a") as log:
                log.write(f"{prefix}\n")
        except Exception as e:
            print(f"-> Failed query for prefix {prefix}: {e}")
            print("Stopping batch execution.")
            break

if __name__ == "__main__":
    main()
