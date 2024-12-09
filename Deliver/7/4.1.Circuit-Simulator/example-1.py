from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram

# Create a simple quantum circuit
qc = QuantumCircuit(2)  # Two qubits
qc.h(0)                 # Apply a Hadamard gate to qubit 0
qc.cx(0, 1)             # Apply a CNOT gate with qubit 0 as control and qubit 1 as target
qc.measure_all()        # Measure all qubits

# Display the circuit
print(qc)

# Initialize the AerSimulator with the tensor_network method
simulator = AerSimulator(device='GPU',method='tensor_network')
#simulator = AerSimulator(device='GPU')

# Transpile the circuit for the simulator
transpiled_qc = transpile(qc, simulator)

# Execute the circuit
job = simulator.run(transpiled_qc, shots=1024)  # 1024 shots for statistical results
result = job.result()

# Get and print the counts (measurement results)
counts = result.get_counts(transpiled_qc)
print("Counts:", counts)

## Visualize the results
#plot_histogram(counts).show()

