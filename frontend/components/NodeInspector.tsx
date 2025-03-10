import React, { useState, useEffect } from 'react';
import { Node } from 'reactflow';
import styles from '../styles/NodeInspector.module.css';

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
        {Object.entries(nodeProperties).map(([key, value]) => (
          <div key={key} className={styles.formGroup}>
            <label htmlFor={`property-${key}`}>{key}:</label>
            <input
              id={`property-${key}`}
              type="text"
              value={value as string}
              onChange={(e) => handlePropertyChange(key, e.target.value)}
            />
          </div>
        ))}
        
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