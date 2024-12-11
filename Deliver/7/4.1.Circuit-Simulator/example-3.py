from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram

class My_Quantum_Circuit_Thing:
    def __init__(self, num_qubits):
        self.num_qubits = num_qubits
        self.circuit = QuantumCircuit(self.num_qubits)
        self.simulator = AerSimulator(device='GPU', method='tensor_network')
        self.transpiled_circuit = None
        self.results = None

    def build_circuit(self, parameters):
        """Builds the quantum circuit.  Now handles more qubits
        and parameters. Uses a simple, repeating pattern for
        demonstration."""

        if len(parameters) < self.num_qubits -1:
             raise ValueError("Insufficient parameters for the circuit.")


        self.circuit.h(0)
        for i in range(self.num_qubits - 1):
            self.circuit.cx(i, i + 1)
            self.circuit.ry(parameters[i], i + 1)  # Using parameters for rotations


        self.circuit.measure_all()


    def run_simulation(self, shots=1024):
        """Transpiles and runs the quantum circuit."""

        self.transpiled_circuit = transpile(self.circuit, self.simulator)
        job = self.simulator.run(self.transpiled_circuit, shots=shots)
        self.results = job.result()

    def get_counts(self):
        """Returns the counts (measurement results)."""

        if self.results is None:
            raise RuntimeError("Circuit has not been run yet.")
        return self.results.get_counts(self.transpiled_circuit)



# Example usage for a 30-qubit circuit:
num_qubits = 30
my_circ = My_Quantum_Circuit_Thing(num_qubits)

# Example parameters – you'll likely want to generate these
# programmatically for a real application.
parameters = [i/100 for i in range(num_qubits-1)]  

my_circ.build_circuit(parameters)
my_circ.run_simulation()
counts = my_circ.get_counts()
print("Counts:", counts)
# plot_histogram(counts).show()  # Can be very large to display


