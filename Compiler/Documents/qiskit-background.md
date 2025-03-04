# Qiskit background: transpile and friends

## Long story short

We first only need to create a `Backend` class. This will be used in
transpilation.

The `Provider` class is for dealing with simulators, authentication etc. Each
actual provider (say IQM, Julich, etc) will write implement a `Provider` class.

As an example of how to create `Backend`, here's example code taken from
[https://docs.quantum.ibm.com/guides/custom-backend](https://docs.quantum.ibm.com/guides/custom-backend)
showing a fake backend and a transpilation snippet, see the folder
`Compiler/examples` and the example below.

Usage:

    % cd Compiler/examples
    % python3 example-transpile.py
    Pre-Transpilation:
    CX gates: 49
    H gates: 50

    ##############################

    Post-Transpilation:
    CZ gates: 142
    ECR gates: 8
    SX gates: 302
    RZ gates: 272

The files:

- `fakeosq_backend.py` contains the `FakeOSQBackend` backend class.

- `example-transpile.py` shows (above) how to use that in actual transpilation.


## Qiskit transpiler

What Qiskit calls transpilation is actually compilation. An elementary example
of how to create, transpile and run a circuit is here

[https://github.com/Qiskit/qiskit-tutorials/blob/master/tutorials/circuits/01_circuit_basics.ipynb](https://github.com/Qiskit/qiskit-tutorials/blob/master/tutorials/circuits/01_circuit_basics.ipynb)

see input `[11]` in the document. The description before the cell is
mistakenly talking about the `execute` function, whereas the cell contains the
transpile-and-run sequence. Also, the 'transpiled' circuit is actually called
`qc_compiled`, which tells us that transpilation is actually compilation.

What is the input and what is output language? The comment in the code says:

> First we have to transpile the quantum circuit to the low-level QASM
> instructions used by the backend.

More detailed explanations of how Qiskit transpilation works are here:

- [Introduction to transpilation](https://docs.quantum.ibm.com/guides/transpile)

- [Transpiler stages](https://docs.quantum.ibm.com/guides/transpiler-stages)

- [Transpiler: qiskit.transpiler](https://docs.quantum.ibm.com/api/qiskit/transpiler#layout-stage)


## We need to implement a (OSQ) backend: IQM example

We can take IQM as an example of how to implement the Python objects needed for
running one's own backend in the Qiskit framework.

- [Qiskit on IQM](https://github.com/iqm-finland/qiskit-on-iqm)

Start with an elementary code example: the following must run

[https://github.com/iqm-finland/qiskit-on-iqm/blob/main/src/iqm/qiskit_iqm/examples/transpile_example.py](https://github.com/iqm-finland/qiskit-on-iqm/blob/main/src/iqm/qiskit_iqm/examples/transpile_example.py)

It's basically this:

    backend = IQMProvider(server_url).get_backend()

    num_qubits = min(backend.num_qubits, 5)  # use at most 5 qubits
    circuit = QuantumCircuit(num_qubits)
    circuit.h(0)
    for i in range(1, num_qubits):
        circuit.cx(0, i)
    circuit.measure_all()

    transpiled_circuit = transpile(circuit, backend)
    counts = backend.run(transpiled_circuit, shots=1000).result().get_counts()

So they have implemented:

- `IQMProvider` in `iqm_provider.py` and `iqm_backend.py`,

which returns a `backend` object. In general, IQM have implemented more, also
transpiler passes,

[https://github.com/iqm-finland/qiskit-on-iqm/tree/main/src/iqm/qiskit_iqm](https://github.com/iqm-finland/qiskit-on-iqm/tree/main/src/iqm/qiskit_iqm)

However, we can start by by implementing the Provider and Backend classes, and
then see if we can add anything into the transpiler.


## Better documentation: IBM on how implement custom Provider, Backend and transpile

IBM has collected all the information on how to write custom transpilers with
custom backends here:

- [Transpile against custom backends](https://docs.quantum.ibm.com/guides/custom-backend)

It also features explicit code examples.


## Chalmers Tergite stack

This is not directly relevant, but let's mention that Chalmers have now open
sourced their stack to run their quantum devices.

- [Tergite](https://tergite.github.io/)

It is mostly the frontend and backend, the frontend to take in jobs from the
user and the backend that sits in front of the QC, which actually runs the jobs
on the device.

