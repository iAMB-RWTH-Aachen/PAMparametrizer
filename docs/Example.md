# Example usage of the PAMParametrizer

## Example 1: Parametrizing the *E. coli* Protein Allocation Model with glucose as sole carbon source
In the first tutorial, we walk over the steps how to set up, run and analyze the result of the PAMparametrizer for the
well studied example of *Escherichia coli*. In case you want to adjust this tutorial to your microbe of study, please
refer to the [PAModelpy documentation](https://iamb-rwth-aachen.github.io/PAModelpy/example) on how to set up the 
PAM for another microbe.

For this entire tutorial, you'll need to load the following packagesL

```python
import os
import pandas as pd
# you want to filter out the warnings to ignore all the times the parametrizer hits infeasibilities, which is pretty annoying
import warnings
warnings.filterwarnings("ignore")

from PAModelpy.configuration import Config

from Modules.PAM_parametrizer import ValidationData, HyperParameters, ParametrizationResults
from Modules.PAM_parametrizer import PAMParametrizer
```

### Step 0: Initiate the parameter set using GotEnzymes
We first need to get our initial parameter set. For this, we use the reactions and proteins which are available in the model,
uniprot and GotEnzymes. In order to parse everything properly, we need to juggle the identifiers of the proteins and reactions.
Please refer to `Scripts/parse_kcat_values_GotEnzymes.ipynb` for an example notebookon how to do this.

### Step 1: Organize the data

### Step 2: Build the Protein Allocation Model

### Step 3: Create the data objects required for the PAMparametrizer
#### i. Get the ValidationData
#### ii. Define the HyperParams
#### iii. ParametrizationResults and FluxResults

### Step 4: Run!

### Step 5: Analyze the Results

## Example 2: Parametrizing the *E. coli* Protein Allocation Model with multiple carbon sources