import time
import random
from datetime import datetime
import os

print("현재 작업 디렉토리:", os.getcwd())
timelap = []
wrongtimes = 0
mistakes = 0
for i in range(10):
    a = random.randint(101, 999)
    b = random.randint(101, 999)
    print(a, "+", b)
    timestart = time.time()
    while True:
        try:
            ans = int(input())
            if ans == a + b:
                elapsed_time = time.time() - timestart
                timelap.append(elapsed_time)
                break
            else:
                wrongtimes += 1
        except ValueError:
            mistakes += 1
            continue
print(timelap)
print(wrongtimes)
print(mistakes)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

with open("timelogs.txt", "a", encoding="utf-8") as f:
    f.write(
        f"{sum(timelap)/len(timelap):.2f}   {wrongtimes}    {mistakes}    {timestamp}\n"
    )
