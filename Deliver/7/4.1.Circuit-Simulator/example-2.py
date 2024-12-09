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
        """
        Builds the quantum circuit based on a list of gate specifications.

        Args:
            parameters: A dummy list argument holding the only parameter
                   that the quantum circuit uses.
        """
        self.circuit.h(0)
        self.circuit.cx(0, 1)
        self.circuit.ry(parameters[0], 1)
        self.circuit.cx(1, 2)
        self.circuit.measure_all()  # Direct measurement of all qubits



    def run_simulation(self, shots=1024):
        """
        Transpiles and runs the quantum circuit.
        """

        self.transpiled_circuit = transpile(self.circuit, self.simulator)
        job = self.simulator.run(self.transpiled_circuit, shots=shots)
        self.results = job.result()

    def get_counts(self):
        """
        Returns the counts (measurement results).
        """

        if self.results is None:
            raise RuntimeError("Circuit has not been run yet.")
        return self.results.get_counts(self.transpiled_circuit)



# Example usage for a more complex circuit:
my_circ = My_Quantum_Circuit_Thing(3)  # 3 qubits

parameters = [0.2345]
my_circ.build_circuit(parameters)  # Constructing circuit


my_circ.run_simulation()
counts = my_circ.get_counts()
print("Counts:", counts)
#plot_histogram(counts).show()


