import subprocess
count = 0
tot=100
for i in range(1, tot):
    process = subprocess.Popen(['python3', 'to_test_reed_soloman.py'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    output = output.decode().strip()
    num1, num2 = map(int, output.split())
    if num1 == num2:
        count += 1

print(f'{count}/{tot-1}')
