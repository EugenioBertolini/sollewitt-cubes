import numpy as np
import itertools
from tqdm import tqdm

sides = 12
minimal_moves = np.array(
    [
        [1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8],
        [5, 9, 6, 1, 0, 8, 10, 2, 4, 11, 7, 3],
        [2, 6, 10, 7, 3, 1, 9, 11, 0, 5, 8, 4],
    ],
    dtype=np.int8,
)
actions = (
    (),
    (0,),
    (0, 0),
    (0, 0, 0),
    (1,),
    (1, 1),
    (1, 1, 1),
    (2,),
    (2, 2),
    (2, 2, 2),
    (0, 2),
    (0, 2, 2),
    (0, 2, 2, 2),
    (0, 0, 2),
    (0, 0, 2, 2, 2),
    (0, 0, 0, 2),
    (0, 0, 0, 2, 2),
    (0, 0, 0, 2, 2, 2),
    (1, 2),
    (1, 2, 2),
    (1, 2, 2, 2),
    (1, 1, 1, 2),
    (1, 1, 1, 2, 2),
    (1, 1, 1, 2, 2, 2),
)

initial_state = np.arange(sides, dtype=np.int8)
all_moves_list = []
for action_sequence in actions:
    current_state = initial_state.copy()
    if not action_sequence == ():
        for move_id in action_sequence:
            current_state = current_state[minimal_moves[move_id]]
    all_moves_list.append(current_state)
all_moves = np.stack(all_moves_list)
print(f"Generated a {all_moves.shape} NumPy array with all move mappings.\n")
print(all_moves)

all_combinations = np.array(list(itertools.product([0, 1], repeat=sides))).astype(
    np.int8
)
idx_done = []
new_cubes = []
for i, this_cube in tqdm(enumerate(all_combinations)):
    if i in idx_done:
        continue

    new_cubes.append(this_cube)
    idx_done.sort()

    for move in all_moves:
        rotated_cube = np.zeros_like(this_cube)
        rotated_cube[move] = this_cube
        for j, other_cube in enumerate(all_combinations[i + 1 :]):
            if all(x == y for x, y in zip(rotated_cube, other_cube)):
                idx_done.append(j + i + 1)

    idx_done = list(set(idx_done))

new_cube_path = "new_cubes.npy"
np.save(new_cube_path, np.array(new_cubes))
