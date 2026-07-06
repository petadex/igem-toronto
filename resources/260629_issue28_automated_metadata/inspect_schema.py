import duckdb

def main():
    print("Initializing DuckDB...")
    con = duckdb.connect()
    
    # Enable S3 access
    con.execute("INSTALL httpfs;")
    con.execute("LOAD httpfs;")
    
    # Configure anonymous access for public S3 bucket
    con.execute("SET s3_region='us-east-1';")
    
    url = "s3://sra-pub-metadata-us-east-1/sra/metadata/20260702_033453_00007_gpbyb_04d49e50-6dd3-4bc0-8636-2f9a03b402da"
    print(f"Reading schema from {url}...")
    try:
        # We can use DESCRIBE to get column names and types
        schema = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{url}') LIMIT 1;").df()
        print("Schema columns:")
        for idx, row in schema.iterrows():
            print(f"  {row['column_name']}: {row['column_type']}")
    except Exception as e:
        print(f"Failed to read schema: {e}")

if __name__ == '__main__':
    main()
