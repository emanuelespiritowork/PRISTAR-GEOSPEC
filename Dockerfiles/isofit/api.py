from fastapi import FastAPI
from pydantic import BaseModel
import subprocess

app = FastAPI()

class Cmd(BaseModel):
    command: str

@app.post("/run")
def run(cmd: Cmd):
    result = subprocess.run(
        cmd.command, shell=True, executable="/bin/bash",
        capture_output=True, text=True
    )
    return {"stdout": result.stdout, "stderr": result.stderr, "returncode": result.returncode}

@app.get("/health")
def health():
    return {"status": "ok"}