from fastapi import FastAPI, UploadFile
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import pandas as pd

app = FastAPI()

# Import R packages
phyloseq = importr('phyloseq')
vegan = importr('vegan')
base = importr('base')

@app.post("/analyze/alpha-diversity")
async def analyze_alpha_diversity(data: UploadFile):
    # Reference the alpha_diversity.Rmd parameters
    # Lines 28-39 show the required parameters
    physeq = import_data(data)
    
    # Run analysis using the R functions
    result = robjects.r('''
        function(physeq) {
            dt <- alpha_diversity(physeq)
            return(dt)
        }
    ''')(physeq)
    
    return {"results": pandas2ri.rpy2py(result)}
