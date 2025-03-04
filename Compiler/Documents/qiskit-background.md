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


## We need to implement a (OSQ) backend



