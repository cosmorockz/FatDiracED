import numpy as np
from header import kinConst, update_intConst

# Define the range of intConst values
intConst_values = np.linspace(-16, 16, 81)

# Loop over intConst values
for new_intConst_value in intConst_values:
    # Update the value of intConst in header.py
    update_intConst(kinConst*new_intConst_value)

    # Run the main code from main.py
    exec(open("main.py").read())

    # Optionally, you can save or process the results here












