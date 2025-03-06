import { useState } from 'react'
import { Box, Button, Card, CardContent, FormControl, InputLabel, MenuItem, Select, Typography } from '@mui/material'

const AVAILABLE_ANALYSES = [
  {
    id: 'alpha_diversity',
    name: 'Alpha Diversity',
    description: 'Measures species diversity within individual samples',
    parameters: {
      indices: ['Richness', 'Shannon', 'Simpson', 'Fisher', 'Pielou']
    }
  },
  {
    id: 'beta_diversity',
    name: 'Beta Diversity', 
    description: 'Measures differences between samples',
    parameters: {
      metrics: ['bray', 'unifrac', 'wunifrac']
    }
  }
]

export const AnalysisForm = () => {
  const [selectedAnalyses, setSelectedAnalyses] = useState<string[]>([])
  const [files, setFiles] = useState<FileList | null>(null)

  const handleSubmit = async () => {
    // Submit to Python API
  }

  return (
    <Card>
      <CardContent>
        <Typography variant="h5">Analysis Configuration</Typography>
        {/* Add form fields */}
      </CardContent>
    </Card>
  )
}
