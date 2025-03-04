from qiskit.visualization import plot_gate_map, plot_coupling_map, plot_circuit_layout
from fakeosq_backend import FakeOSQBackend

backend = FakeOSQBackend(3, 3)

target = backend.target
coupling_map_backend = target.build_coupling_map()

coordinates = [
    (3, 1),
    (3, -1),
    (2, -2),
    (1, 1),
    (0, 0),
    (-1, -1),
    (-2, 2),
    (-3, 1),
    (-3, -1),
    (2, 1),
    (1, -1),
    (-1, 1),
    (-2, -1),
    (3, 0),
    (2, -1),
    (0, 1),
    (0, -1),
    (-2, 1),
    (-3, 0),
]

single_qubit_coordinates = []
total_qubit_coordinates = []


for coordinate in coordinates:
    total_qubit_coordinates.append(coordinate)

for coordinate in coordinates:
    total_qubit_coordinates.append(
        (-1 * coordinate[0] + 1, coordinate[1] + 4)
    )

for coordinate in coordinates:
    total_qubit_coordinates.append((coordinate[0], coordinate[1] + 8))


line_colors = ["#adaaab" for edge in coupling_map_backend.get_edges()]
ecr_edges = []

# Get tuples for the edges which have an ecr instruction attached
for instruction in target.instructions:
    if instruction[0].name == "ecr":
        ecr_edges.append(instruction[1])

for i, edge in enumerate(coupling_map_backend.get_edges()):
    if edge in ecr_edges:
        line_colors[i] = "#000000"

print(backend.name)
plot_gate_map(
    backend,
    plot_directed=True,
    qubit_coordinates=total_qubit_coordinates,
    line_color=line_colors,
)
