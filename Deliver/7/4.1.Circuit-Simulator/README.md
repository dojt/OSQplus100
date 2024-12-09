# DOT's folder for collecting the output for Task 7.4.1


### Step 1: Set up your Python environment

Only once, do this:

```bash
python3 -m venv name_of_the_venv # `name_of_the_venv` better be short!
```

Everytime, do this:

```bash
source name_of_the_venv/bin/activate
```

Now, do your thing (steps 2,3,...).  When you're done, every time, do this:

```bash
deactivate
```

### Step 2: Install Qiskit etc

This step you only have to do once.

```bash
pip install qiskit
pip install qiskit-aer
pip install qiskit-aer-gpu
```

Run "example-1.py" to see if it works.  If it doesn't, move to Réunion and train to be a donkey drover.
Continue with step 3.

### Step 3: Run your code, using `AerSimulator(device="GPU",method="tensor_network")`

That's basically it.  But the examples files have more ... examples.
