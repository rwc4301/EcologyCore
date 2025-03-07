import { useCallback } from 'react';
import ReactFlow, {
  Background,
  Controls,
  MiniMap,
  ReactFlowProvider,
  Panel
} from 'reactflow';
import 'reactflow/dist/style.css';

const initialNodes = [
  {
    id: '1',
    type: 'input',
    data: { label: 'Input Node' },
    position: { x: 250, y: 25 },
  }
];

const initialEdges = [];

export default function WorkflowEditor() {
  const onNodesChange = useCallback((changes: any) => {
    // Will implement node changes later
    console.log('Node changed:', changes);
  }, []);

  const onEdgesChange = useCallback((changes: any) => {
    // Will implement edge changes later
    console.log('Edge changed:', changes);
  }, []);

  const onConnect = useCallback((connection: any) => {
    // Will implement connection handling later
    console.log('New connection:', connection);
  }, []);

  return (
    <div style={{ width: '100%', height: '80vh' }}>
      <ReactFlowProvider>
        <ReactFlow
          nodes={initialNodes}
          edges={initialEdges}
          onNodesChange={onNodesChange}
          onEdgesChange={onEdgesChange}
          onConnect={onConnect}
          fitView
        >
          <Background />
          <Controls />
          <MiniMap />
          <Panel position="top-left">
            <h3>Workflow Editor</h3>
          </Panel>
        </ReactFlow>
      </ReactFlowProvider>
    </div>
  );
} 