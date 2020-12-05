from global_relaxation import global_relaxation
from local_relaxation import local_relaxation
import datetime

if __name__ == "__main__":
    start = datetime.datetime.now()
    global_relaxation()
    local_relaxation()
    end = datetime.datetime.now()
    print(f'time: {end - start}')