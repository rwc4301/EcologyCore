import { memo } from 'react';
import { Handle, Position, NodeProps } from 'reactflow';
import styles from '../styles/CustomNodes.module.css';

// Define the data structure for MultiPortNode
interface MultiPortNodeData {
  label: string;
  inputs?: Array<{
    id: string;
    label: string;
    type?: string;
    position?: number; // 0-100 percentage position
  }>;
  outputs?: Array<{
    id: string;
    label: string;
    type?: string;
    position?: number; // 0-100 percentage position
  }>;
}

/**
 * A custom node component that supports multiple inputs and outputs
 */
export const MultiPortNode = memo(({ data, id }: NodeProps<MultiPortNodeData>) => {
  const nodeData = data || { label: 'Node' };
  const inputs = nodeData.inputs || [];
  const outputs = nodeData.outputs || [];
  
  // Convert position percentage to CSS positioning
  const getHandlePosition = (position: number | undefined, index: number, total: number) => {
    if (position !== undefined) {
      // Use the provided position (as a percentage)
      return `${Math.min(Math.max(position, 0), 100)}%`;
    } else {
      // Distribute evenly
      const step = 100 / (total + 1);
      return `${step * (index + 1)}%`;
    }
  };
  
  return (
    <div className={styles.multiPortNode}>
      <div className={styles.nodeHeader}>{nodeData.label}</div>
      
      {/* Input Handles */}
      {inputs.map((input, index) => (
        <div 
          key={`input-${input.id}`} 
          className={styles.portLabel}
          style={{ 
            position: 'absolute', 
            left: '0px',
            top: getHandlePosition(input.position, index, inputs.length),
            transform: 'translateY(-50%)'
          }}
        >
          <Handle
            type="target"
            position={Position.Left}
            id={input.id}
            className={styles.handle}
            style={{ left: '-8px' }}
          />
          <span className={styles.portName}>
            {input.label}
            {input.type && <span className={styles.portType}>{input.type}</span>}
          </span>
        </div>
      ))}
      
      {/* Output Handles */}
      {outputs.map((output, index) => (
        <div 
          key={`output-${output.id}`} 
          className={styles.portLabel}
          style={{ 
            position: 'absolute', 
            right: '0px',
            top: getHandlePosition(output.position, index, outputs.length),
            transform: 'translateY(-50%)'
          }}
        >
          <span className={styles.portName}>
            {output.label}
            {output.type && <span className={styles.portType}>{output.type}</span>}
          </span>
          <Handle
            type="source"
            position={Position.Right}
            id={output.id}
            className={styles.handle}
            style={{ right: '-8px' }}
          />
        </div>
      ))}
      
      {/* If no inputs or outputs are defined, show default handles for compatibility */}
      {inputs.length === 0 && (
        <Handle type="target" position={Position.Left} className={styles.defaultHandle} />
      )}
      {outputs.length === 0 && (
        <Handle type="source" position={Position.Right} className={styles.defaultHandle} />
      )}
    </div>
  );
});

// Define node types object to use with ReactFlow
export const nodeTypes = {
  multiPort: MultiPortNode,
};