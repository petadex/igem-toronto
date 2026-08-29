

import boto3
import textwrap

REGION = "us-east-1"
AMI_ID = "ami-XXXXXXXXXXXXXXXXX"          # your prebuilt Sitetack AMI
INSTANCE_TYPE = "g6.xlarge"
KEY_NAME = "your-keypair-name"             # for SSH fallback/debugging
SECURITY_GROUP_IDS = ["sg-XXXXXXXXXXXXXXXXX"]
SUBNET_ID = "subnet-XXXXXXXXXXXXXXXXX"
IAM_INSTANCE_PROFILE = "sitetack-batch-role"   # grants S3 read/write

S3_INPUT_PREFIX = "s3://your-bucket/petadex/sequences/"
S3_MODEL_PREFIX = "s3://your-bucket/petadex/models/"
S3_OUTPUT_PREFIX = "s3://your-bucket/petadex/predictions/"

MODEL_NAMES = [
    "model_01", "model_02", "model_03", "model_04", "model_05",
    "model_06", "model_07", "model_08", "model_09", "model_10",
    "model_11", "model_12", "model_13",
]

# Max Spot price you're willing to pay per hour (safety cap).
MAX_SPOT_PRICE = "1"

# ------------------------------------------------------------------------


def build_user_data(model_name: str) -> str:
    """
    Generates the per-instance startup script. It:
      1. Pulls that model's .h5 file and the shared input sequences from S3
      2. Activates the sitetack conda env
      3. Runs the batch prediction for this model only
      4. Uploads results back to S3
      5. Shuts the instance down (terminates, since we set
         InstanceInitiatedShutdownBehavior=terminate at launch)
    Logs go to /var/log/sitetack-run.log on the instance for debugging
    before it terminates (check CloudWatch Logs if you wire that up,
    since the instance disappears after shutdown).
    """
    script = f"""\
        #!/bin/bash
        set -euo pipefail
        exec > /var/log/sitetack-run.log 2>&1

        echo "=== Starting Sitetack batch for {model_name} at $(date) ==="

        WORKDIR=/home/ubuntu/sitetack-run
        mkdir -p "$WORKDIR"
        cd "$WORKDIR"

        # Pull this model's weights and the shared input sequences
        aws s3 cp "{S3_MODEL_PREFIX}{model_name}.h5" ./{model_name}.h5
        aws s3 sync "{S3_INPUT_PREFIX}" ./sequences/

        # Activate conda env (adjust path if conda isn't in this location on your AMI)
        source /home/ubuntu/miniconda3/etc/profile.d/conda.sh
        conda activate sitetack

        # Run the batch prediction for this model
        # Replace this with however you actually invoke Sitetack for one model
        python run_sitetack_batch.py \\
            --model ./{model_name}.h5 \\
            --input ./sequences/ \\
            --output ./results_{model_name}.csv

        # Push results back to S3
        aws s3 cp ./results_{model_name}.csv "{S3_OUTPUT_PREFIX}results_{model_name}.csv"

        echo "=== Finished {model_name} at $(date), shutting down ==="
        shutdown -h now
    """
    return textwrap.dedent(script)


def launch_spot_instance(ec2_client, model_name: str):
    user_data = build_user_data(model_name)

    response = ec2_client.run_instances(
        ImageId=AMI_ID,
        InstanceType=INSTANCE_TYPE,
        KeyName=KEY_NAME,
        MinCount=1,
        MaxCount=1,
        SecurityGroupIds=SECURITY_GROUP_IDS,
        SubnetId=SUBNET_ID,
        IamInstanceProfile={"Name": IAM_INSTANCE_PROFILE},
        UserData=user_data,
        InstanceInitiatedShutdownBehavior="terminate",
        InstanceMarketOptions={
            "MarketType": "spot",
            "SpotOptions": {
                "MaxPrice": MAX_SPOT_PRICE,
                "SpotInstanceType": "one-time",
                "InstanceInterruptionBehavior": "terminate",
            },
        },
        TagSpecifications=[
            {
                "ResourceType": "instance",
                "Tags": [
                    {"Key": "Name", "Value": f"sitetack-{model_name}"},
                    {"Key": "Project", "Value": "PETadex"},
                    {"Key": "Model", "Value": model_name},
                ],
            }
        ],
    )
    instance_id = response["Instances"][0]["InstanceId"]
    print(f"Launched {model_name} -> {instance_id}")
    return instance_id


def main():
    ec2 = boto3.client("ec2", region_name=REGION)
    launched = {}
    for model_name in MODEL_NAMES:
        instance_id = launch_spot_instance(ec2, model_name)
        launched[model_name] = instance_id

    print("\nAll instances launched:")
    for model_name, instance_id in launched.items():
        print(f"  {model_name}: {instance_id}")

    print(
        "\nEach instance will self-terminate after uploading its results to "
        f"{S3_OUTPUT_PREFIX}. Poll with:\n"
        f"  aws s3 ls {S3_OUTPUT_PREFIX}\n"
        "or check instance state with:\n"
        "  aws ec2 describe-instances --filters "
        '"Name=tag:Project,Values=PETadex" '
        '"Name=instance-state-name,Values=running,pending"'
    )


if __name__ == "__main__":
    main()

