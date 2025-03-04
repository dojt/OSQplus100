# Qiskit background: transpile and friends

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


## Chalmers Tergite stack

This is not directly relevant, but let's mention that Chalmers have now open
sourced their stack to run their quantum devices.

- [Tergite](https://tergite.github.io/)

It is mostly the frontend and backend, the frontend to take in jobs from the
user and the backend that sits in front of the QC, which actually runs the jobs
on the device.

