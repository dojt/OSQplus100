from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
import random

class My_Quantum_Circuit_Thing:
    def __init__(self, num_qubits):
        self.num_qubits = num_qubits
        self.circuit = QuantumCircuit(self.num_qubits)
        self.simulator = AerSimulator(device='GPU', method='tensor_network')
        self.transpiled_circuit = None
        self.results = None

    def build_circuit(self, num_gates):
        """Builds the quantum circuit with at least num_gates gates."""
        
        if num_gates < 10 * self.num_qubits:
            raise ValueError(f"Number of gates must be at least 10 times the number of qubits ({10 * self.num_qubits}).")

        # Apply an initial layer of Hadamard gates
        self.circuit.h(range(self.num_qubits))

        gates_applied = 0
        while gates_applied < num_gates:
            for i in range(self.num_qubits):
                # Randomly choose a single-qubit gate (H, RX, RY, RZ)
                gate_choice = random.choice(['h', 'rx', 'ry', 'rz'])
                if gate_choice == 'h':
                    self.circuit.h(i)
                elif gate_choice == 'rx':
                    self.circuit.rx(random.uniform(0, 2 * 3.14159), i)  # Random rotation angle
                elif gate_choice == 'ry':
                    self.circuit.ry(random.uniform(0, 2 * 3.14159), i)
                elif gate_choice == 'rz': # Or just else: would work too
                    self.circuit.rz(random.uniform(0, 2 * 3.14159), i) 
                gates_applied += 1

                if gates_applied >= num_gates:
                    break  #Exit if enough gates have been applied

                # Add a two-qubit gate if there are at least 2 qubits and we haven't reached the target
                if self.num_qubits >=2 and gates_applied < num_gates:
                    control_qubit = random.randint(0, self.num_qubits - 1)
                    target_qubit = random.randint(0, self.num_qubits - 1)
                    while target_qubit == control_qubit: # Ensure different qubits
                        target_qubit = random.randint(0, self.num_qubits-1)
                    self.circuit.cx(control_qubit, target_qubit)
                    gates_applied += 1

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



# Example usage:
num_qubits = 30
num_gates = 10 * num_qubits   # Set number of gates to be added

my_circ = My_Quantum_Circuit_Thing(num_qubits)
my_circ.build_circuit(num_gates)
my_circ.run_simulation()
counts = my_circ.get_counts()
print("Counts:", counts)

