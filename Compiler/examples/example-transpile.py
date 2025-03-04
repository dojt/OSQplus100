from fakeosq_backend import FakeOSQBackend
from qiskit import QuantumCircuit, transpile
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

backend = FakeOSQBackend(3, 3)

num_qubits = 50
ghz = QuantumCircuit(num_qubits)
ghz.h(range(num_qubits))
ghz.cx(0, range(1, num_qubits))
op_counts = ghz.count_ops()

print("Pre-Transpilation: ")
print(f"CX gates: {op_counts['cx']}")
print(f"H gates: {op_counts['h']}")
print()
print(30 * "#", "\n")

pm = generate_preset_pass_manager(optimization_level=3, backend=backend)
transpiled_ghz = pm.run(ghz)
op_counts = transpiled_ghz.count_ops()

print("Post-Transpilation: ")
print(f"CZ gates: {op_counts['cz']}")
print(f"ECR gates: {op_counts['ecr']}")
print(f"SX gates: {op_counts['sx']}")
print(f"RZ gates: {op_counts['rz']}")
