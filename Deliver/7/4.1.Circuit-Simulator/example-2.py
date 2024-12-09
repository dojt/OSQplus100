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
            gates: A list of tuples, where each tuple represents a gate.
                   The tuple format is (gate_name, target_qubits, params).
                   For example: ('h', [0], []), ('cx', [0, 1], []), ('rx', [1], [0.5]).
        """

        gates = [
           ('h', [0], []),  # Hadamard on qubit 0
            ('cx', [0, 1], []),  # CNOT on qubits 0 and 1
            ('ry', [1], [parameters[0]]),  # Rotation around y-axis on qubit 1
            ('cx', [1,2], []), # CNOT on qubits 1 and 2
            ('measure', [0,1,2], [])
        ]

        for gate_name, target_qubits, params in gates:

            if hasattr(self.circuit,gate_name):
                gate_method = getattr(self.circuit, gate_name)
                gate_method(*params, *target_qubits)  # Apply gate

            else:
                raise ValueError(f"Gate {gate_name} is not supported.")



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


