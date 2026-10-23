import os
import time

# Directory paths
input_dir = "/home/belghith1/Exhaustive_Z4/Cyclic/Input"
output_dir = "/home/belghith1/Exhaustive_Z4/Cyclic/Output"
utilities_dir = "/home/belghith1/Exhaustive_Z4/Cyclic/Utilities"
template_path = "/home/belghith1/Exhaustive_Z4/Cyclic/Scripts/template.m"

# Read the template file
with open(template_path, "r") as file:
    magma_script_template = file.read()

def generate_magma_scripts():
    # Create the input directory if it does not exist
    os.makedirs(input_dir, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(utilities_dir, exist_ok=True)

    # Create Magma scripts for odd lengths from 5 to 97
    for n in range(5, 98, 2):
        print(f"Generating script for length {n}")
        script_content = magma_script_template.format(n=n, output_dir=output_dir)
        script_path = os.path.join(input_dir, f"Len{n}.m")
        with open(script_path, "w") as script_file:
            script_file.write(script_content)
        print(f"Script for length {n} saved to {script_path}")

    # Utility script to run all Magma scripts
    utility_script_path = os.path.join(utilities_dir, "run_all.sh")
    tracking_file_path = os.path.join(utilities_dir, "process_tracking.txt")

    print(f"Creating utility script at {utility_script_path}")

    with open(utility_script_path, "w") as utility_script:
        utility_script.write("#!/bin/bash\n\n")
        utility_script.write(f"echo 'Process Name,PID,Start Time,End Time,Duration,Memory Usage (MB)' > {tracking_file_path}\n\n")
        for n in range(5, 98, 2):
            script_path = os.path.join(input_dir, f"Len{n}.m")
            log_path = os.path.join(output_dir, f"Len{n}_log.txt")
            utility_script.write(f"""
start_time=$(date +%s)
echo "Running Len{n}.m"
magma {script_path} > {log_path} 2>&1 &
pid=$!
echo "Started process $pid for Len{n}.m"
while kill -0 $pid 2>/dev/null; do
    sleep 1
done
end_time=$(date +%s)
duration=$((end_time - start_time))
mem_usage=$(ps -o rss= -p $pid 2>/dev/null)
mem_usage=${{mem_usage:-0}} # Default to 0 if mem_usage is empty
echo "Len{n}.m,$pid,$start_time,$end_time,$duration,$mem_usage" >> {tracking_file_path}
""")
        print(f"Added script for length {n} to utility script")

    # Make the utility script executable
    os.chmod(utility_script_path, 0o755)
    print(f"Utility script saved to {utility_script_path}")

if __name__ == "__main__":
    generate_magma_scripts()
    print("All scripts have been generated successfully.")
