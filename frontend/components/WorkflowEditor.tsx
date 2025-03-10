import { useCallback, useState } from 'react';
import ReactFlow, {
  Background,
  Controls,
  MiniMap,
  ReactFlowProvider,
  Panel,
  Node,
  useNodesState,
  useEdgesState,
  NodeChange,
  EdgeChange,
  Connection,
  addEdge
} from 'reactflow';
import 'reactflow/dist/style.css';
import NodeInspector from './NodeInspector';
import { nodeTypes } from './CustomNodes';

const initialNodes: Node[] = [
  {
    id: '1',
    type: 'input',
    data: { label: 'Input Node' },
    position: { x: 250, y: 25 },
  },
  {
    id: '2',
    type: 'default',
    data: { label: 'Process Node', description: 'Processes data from input' },
    position: { x: 250, y: 125 },
  },
  {
    id: '3',
    type: 'output',
    data: { label: 'Output Node', output: 'final result' },
    position: { x: 250, y: 225 },
  },
  {
    id: '4',
    type: 'multiPort',
    data: { 
      label: 'Multi-Port Node',
      inputs: [
        { id: 'input1', label: 'Input 1', type: 'number' },
        { id: 'input2', label: 'Input 2', type: 'string' }
      ],
      outputs: [
        { id: 'output1', label: 'Output 1', type: 'array' },
        { id: 'output2', label: 'Output 2', type: 'object' }
      ]
    },
    position: { x: 500, y: 125 },
  }
];

const initialEdges = [];

export default function WorkflowEditor() {
  // Use the React Flow state hooks
  const [nodes, setNodes, onNodesChange] = useNodesState(initialNodes);
  const [edges, setEdges, onEdgesChange] = useEdgesState(initialEdges);
  const [selectedNode, setSelectedNode] = useState<Node | null>(null);

  // Handle node selection
  const onNodeClick = useCallback((event: React.MouseEvent, node: Node) => {
    setSelectedNode(node);
  }, []);

  // Handle background click to deselect node
  const onPaneClick = useCallback(() => {
    setSelectedNode(null);
  }, []);

  // Handle node updates from the inspector
  const onNodeUpdate = useCallback((nodeId: string, newData: any) => {
    setNodes((nds) =>
      nds.map((node) => {
        if (node.id === nodeId) {
          // Update the node with new data
          const updatedNode = {
            ...node,
            data: {
              ...newData
            },
          };
          
          // Also update the selected node reference
          setSelectedNode(updatedNode);
          
          return updatedNode;
        }
        return node;
      })
    );
  }, [setNodes]);

  // Handle new connections
  const onConnect = useCallback(
    (params: Connection) => setEdges((eds) => addEdge(params, eds)),
    [setEdges]
  );

  return (
    <div style={{ width: '100%', height: '80vh', position: 'relative' }}>
      <ReactFlowProvider>
        <ReactFlow
          nodes={nodes}
          edges={edges}
          nodeTypes={nodeTypes}
          onNodesChange={onNodesChange}
          onEdgesChange={onEdgesChange}
          onConnect={onConnect}
          onNodeClick={onNodeClick}
          onPaneClick={onPaneClick}
          fitView
        >
          <Background />
          <Controls />
          <MiniMap />
          <Panel position="top-left">
            <h3>Workflow Editor</h3>
          </Panel>
          <Panel position="top-right" style={{ margin: '10px' }}>
            <NodeInspector 
              selectedNode={selectedNode} 
              onNodeUpdate={onNodeUpdate} 
            />
          </Panel>
        </ReactFlow>
      </ReactFlowProvider>
    </div>
  );
}