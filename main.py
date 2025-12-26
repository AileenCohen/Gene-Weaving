import os
import subprocess

def launch():
    # Points to the app.py location
    path = os.path.join("src", "gui", "app.py")
    subprocess.run(["streamlit", "run", path])

if __name__ == "__main__":
    launch()