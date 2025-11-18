import subprocess

# Loop through the numbers 21 to 50
for num in range(572, 646):
    # Define the command to submit the job
    command = f"sbatch 4_{num}_docking.sh"
    
    # Execute the command
    subprocess.run(command, shell=True)
    
    print(f"Submitted job 4_{num}_docking.sh")