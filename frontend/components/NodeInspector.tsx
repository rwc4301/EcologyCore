import React, { useState, useEffect } from 'react';
import { Node } from 'reactflow';
import styles from '../styles/NodeInspector.module.css';
import nodeStyles from '../styles/CustomNodes.module.css';

interface NodeInspectorProps {
  selectedNode: Node | null;
  onNodeUpdate: (nodeId: string, newData: any) => void;
}

export default function NodeInspector({ selectedNode, onNodeUpdate }: NodeInspectorProps) {
  const [nodeProperties, setNodeProperties] = useState<Record<string, any>>({});
  
  useEffect(() => {
    if (selectedNode) {
      // Extract node data properties for editing
      setNodeProperties(selectedNode.data || {});
    } else {
      setNodeProperties({});
    }
  }, [selectedNode]);

  const handlePropertyChange = (key: string, value: any) => {
    const updatedProperties = {
      ...nodeProperties,
      [key]: value
    };
    setNodeProperties(updatedProperties);
  };

  const handleSave = () => {
    if (selectedNode) {
      onNodeUpdate(selectedNode.id, { ...selectedNode.data, ...nodeProperties });
    }
  };
  
  const isMultiPortNode = selectedNode?.type === 'multiPort';

  // Handle input port operations
  const addInputPort = () => {
    const inputs = [...(nodeProperties.inputs || [])];
    const newId = `input_${inputs.length + 1}`;
    inputs.push({
      id: newId,
      label: `Input ${inputs.length + 1}`,
      type: 'any'
    });
    handlePropertyChange('inputs', inputs);
  };

  const updateInputPort = (index: number, field: string, value: any) => {
    const inputs = [...(nodeProperties.inputs || [])];
    inputs[index] = { ...inputs[index], [field]: value };
    handlePropertyChange('inputs', inputs);
  };

  const removeInputPort = (index: number) => {
    const inputs = [...(nodeProperties.inputs || [])];
    inputs.splice(index, 1);
    handlePropertyChange('inputs', inputs);
  };

  // Handle output port operations
  const addOutputPort = () => {
    const outputs = [...(nodeProperties.outputs || [])];
    const newId = `output_${outputs.length + 1}`;
    outputs.push({
      id: newId,
      label: `Output ${outputs.length + 1}`,
      type: 'any'
    });
    handlePropertyChange('outputs', outputs);
  };

  const updateOutputPort = (index: number, field: string, value: any) => {
    const outputs = [...(nodeProperties.outputs || [])];
    outputs[index] = { ...outputs[index], [field]: value };
    handlePropertyChange('outputs', outputs);
  };

  const removeOutputPort = (index: number) => {
    const outputs = [...(nodeProperties.outputs || [])];
    outputs.splice(index, 1);
    handlePropertyChange('outputs', outputs);
  };

  if (!selectedNode) {
    return (
      <div className={styles.inspector}>
        <h3>Node Inspector</h3>
        <p>Select a node to view and edit its properties</p>
      </div>
    );
  }

  return (
    <div className={styles.inspector}>
      <h3>Node Inspector</h3>
      <div className={styles.nodeInfo}>
        <p><strong>ID:</strong> {selectedNode.id}</p>
        <p><strong>Type:</strong> {selectedNode.type}</p>
      </div>
      
      <h4>Properties</h4>
      <div className={styles.propertiesForm}>
        {/* Filter out inputs/outputs from basic properties if it's a multiport node */}
        {Object.entries(nodeProperties)
          .filter(([key]) => !isMultiPortNode || (key !== 'inputs' && key !== 'outputs'))
          .map(([key, value]) => {
            // Skip complex objects/arrays in simple view
            if (typeof value === 'object' && value !== null) {
              return null;
            }
            
            return (
              <div key={key} className={styles.formGroup}>
                <label htmlFor={`property-${key}`}>{key}:</label>
                <input
                  id={`property-${key}`}
                  type="text"
                  value={value as string}
                  onChange={(e) => handlePropertyChange(key, e.target.value)}
                />
              </div>
            );
          })}
        
        {/* Add new property section */}
        <details className={styles.addProperty}>
          <summary>Add New Property</summary>
          <AddNewProperty 
            onAdd={(key, value) => {
              if (key.trim()) {
                handlePropertyChange(key, value);
              }
            }} 
          />
        </details>
      </div>
      
      {/* Multi-port node specific controls */}
      {isMultiPortNode && (
        <div className={nodeStyles.portControls}>
          {/* Input Ports Section */}
          <div className={nodeStyles.portSection}>
            <div className={nodeStyles.portHeading}>
              <h4>Input Ports</h4>
              <button onClick={addInputPort}>Add Input</button>
            </div>
            
            {nodeProperties.inputs?.map((input: any, index: number) => (
              <div key={`input-${index}`} className={nodeStyles.portItem}>
                <button 
                  className={nodeStyles.removeButton}
                  onClick={() => removeInputPort(index)}
                >
                  ×
                </button>
                
                <div className={nodeStyles.portForm}>
                  <div className={nodeStyles.formGroup}>
                    <label>ID:</label>
                    <input
                      type="text"
                      value={input.id}
                      onChange={(e) => updateInputPort(index, 'id', e.target.value)}
                    />
                  </div>
                  
                  <div className={nodeStyles.formGroup}>
                    <label>Label:</label>
                    <input
                      type="text"
                      value={input.label}
                      onChange={(e) => updateInputPort(index, 'label', e.target.value)}
                    />
                  </div>
                  
                  <div className={nodeStyles.formGroup}>
                    <label>Type:</label>
                    <input
                      type="text"
                      value={input.type || ''}
                      onChange={(e) => updateInputPort(index, 'type', e.target.value)}
                    />
                  </div>
                  
                  <div className={nodeStyles.formGroup}>
                    <label>Position (0-100%):</label>
                    <input
                      type="number"
                      min="0"
                      max="100"
                      value={input.position || ''}
                      onChange={(e) => updateInputPort(index, 'position', parseInt(e.target.value))}
                    />
                  </div>
                </div>
              </div>
            ))}
            
            {(!nodeProperties.inputs || nodeProperties.inputs.length === 0) && (
              <p>No input ports defined. Add one to connect to this node.</p>
            )}
          </div>
          
          {/* Output Ports Section */}
          <div className={nodeStyles.portSection}>
            <div className={nodeStyles.portHeading}>
              <h4>Output Ports</h4>
              <button onClick={addOutputPort}>Add Output</button>
            </div>
            
            {nodeProperties.outputs?.map((output: any, index: number) => (
              <div key={`output-${index}`} className={nodeStyles.portItem}>
                <button 
                  className={nodeStyles.removeButton}
                  onClick={() => removeOutputPort(index)}
                >
                  ×
                </button>
                
                <div className={nodeStyles.portForm}>
                  <div className={nodeStyles.formGroup}>
                    <label>ID:</label>
                    <input
                      type="text"
                      value={output.id}
                      onChange={(e) => updateOutputPort(index, 'id', e.target.value)}
                    />
                  </div>
                  
                  <div className={nodeStyles.formGroup}>
                    <label>Label:</label>
                    <input
                      type="text"
                      value={output.label}
                      onChange={(e) => updateOutputPort(index, 'label', e.target.value)}
                    />
                  </div>
                  
                  <div className={nodeStyles.formGroup}>
                    <label>Type:</label>
                    <input
                      type="text"
                      value={output.type || ''}
                      onChange={(e) => updateOutputPort(index, 'type', e.target.value)}
                    />
                  </div>
                  
                  <div className={nodeStyles.formGroup}>
                    <label>Position (0-100%):</label>
                    <input
                      type="number"
                      min="0"
                      max="100"
                      value={output.position || ''}
                      onChange={(e) => updateOutputPort(index, 'position', parseInt(e.target.value))}
                    />
                  </div>
                </div>
              </div>
            ))}
            
            {(!nodeProperties.outputs || nodeProperties.outputs.length === 0) && (
              <p>No output ports defined. Add one to connect from this node.</p>
            )}
          </div>
        </div>
      )}
      
      <button 
        className={styles.saveButton}
        onClick={handleSave}
      >
        Save Changes
      </button>
    </div>
  );
}

function AddNewProperty({ onAdd }: { onAdd: (key: string, value: any) => void }) {
  const [newKey, setNewKey] = useState('');
  const [newValue, setNewValue] = useState('');
  
  const handleAdd = () => {
    if (newKey.trim()) {
      onAdd(newKey, newValue);
      setNewKey('');
      setNewValue('');
    }
  };
  
  return (
    <div className={styles.addPropertyForm}>
      <div className={styles.formGroup}>
        <label htmlFor="new-property-key">Property Name:</label>
        <input
          id="new-property-key"
          type="text"
          value={newKey}
          onChange={(e) => setNewKey(e.target.value)}
          placeholder="Enter property name"
        />
      </div>
      <div className={styles.formGroup}>
        <label htmlFor="new-property-value">Property Value:</label>
        <input
          id="new-property-value"
          type="text"
          value={newValue}
          onChange={(e) => setNewValue(e.target.value)}
          placeholder="Enter property value"
        />
      </div>
      <button onClick={handleAdd}>Add Property</button>
    </div>
  );
}