# Push vs Poll: Real-time UI Updates Architecture

## Current State

The frontend uses adaptive polling to keep the UI synchronized with backend state:

| Context | Active Interval | Idle Interval |
|---------|-----------------|---------------|
| Job tree | 3s | 30s |
| Report/Directory/Logs | 5s | 0 (disabled) |

A **grace period mechanism** handles the race condition where a job is submitted but the database status hasn't updated yet, ensuring fast polling continues for 20 seconds after job submission.

## The Push Alternative

### Electron (IPC)

In the Electron app, Inter-Process Communication could push state changes directly to the renderer:

```
Backend (Django)
    → Signal emitted on job state change / file registration
    → Main process receives via bridge
    → IPC message to renderer
    → SWR mutate() triggers refresh
```

**Advantages:**
- Near-instant UI updates
- Reduced network traffic
- Lower CPU usage when idle

### Web Deployment

For web, equivalent push mechanisms:
- **WebSockets** - Full duplex, Django Channels
- **Server-Sent Events (SSE)** - Simpler, one-way push
- **Polling fallback** - Current implementation serves as baseline

## Signal/Slot as Foundation

The existing `HierarchicalObject` Signal/Slot mechanism (`baselayer/`) provides a natural integration point:

```python
# In job execution code
job.status = JobStatus.RUNNING
job.statusChanged.emit()  # Signal already exists

# New: Bridge to transport layer
@statusChanged.connect
def notify_frontend(job):
    if is_electron():
        ipc_bridge.send('job:status', job.id, job.status)
    elif has_websocket():
        ws_channel.send_json({'type': 'job:status', ...})
```

The signals are **transport-agnostic** - the same signal can trigger IPC in Electron or be queued for WebSocket/SSE in web mode.

## Hybrid Architecture

A practical implementation would layer push on top of polling:

```
┌─────────────────────────────────────────────────┐
│                   Frontend                       │
├─────────────────────────────────────────────────┤
│  SWR Cache                                       │
│    ↑                                             │
│    ├── Push: mutate() on IPC/WS message         │
│    ├── Poll: Adaptive interval (fallback)        │
│    └── Grace period: Recently started jobs       │
└─────────────────────────────────────────────────┘
                        ↑
┌─────────────────────────────────────────────────┐
│              Transport Layer                     │
├─────────────────────────────────────────────────┤
│  Electron: IPC via main process                  │
│  Web: WebSocket / SSE / Poll                     │
└─────────────────────────────────────────────────┘
                        ↑
┌─────────────────────────────────────────────────┐
│                   Backend                        │
├─────────────────────────────────────────────────┤
│  Django + Signal/Slot (HierarchicalObject)       │
│    - Job status changes                          │
│    - File registrations                          │
│    - Subjob creation                             │
└─────────────────────────────────────────────────┘
```

## Implementation Considerations

### Events to Push

| Event | Trigger | Payload |
|-------|---------|---------|
| Job status change | `job.save()` when status changes | `{jobId, status, projectId}` |
| File registered | `File.save()` | `{fileId, jobId, projectId}` |
| Subjob created | Child job creation | `{jobId, parentId, projectId}` |
| Validation updated | Parameter change | `{jobId}` |

### Frontend Handling

```typescript
// Hypothetical IPC listener in Electron
useEffect(() => {
  const unsubscribe = ipcRenderer.on('job:status', (jobId, status) => {
    // Trigger SWR revalidation for affected endpoints
    mutateJobTree();
    if (currentJobId === jobId) {
      mutateJob();
      mutateReport();
    }
  });
  return unsubscribe;
}, []);
```

### Graceful Degradation

The polling infrastructure provides automatic fallback:
- Push connection lost → polling continues at adaptive intervals
- Grace period mechanism handles edge cases
- No user-visible degradation, just slightly delayed updates

## Trade-offs

| Aspect | Push | Poll |
|--------|------|------|
| Latency | ~Instant | 3-30s |
| Complexity | Higher (connection state, reconnection) | Lower |
| Server load | Lower (event-driven) | Higher (regular requests) |
| Debugging | Harder (stateful) | Easier (stateless requests) |
| Web compatibility | Requires WebSocket/SSE | Works everywhere |

## Recommendation

1. **Phase 1 (Current):** Polling with adaptive intervals and grace period - provides a solid, debuggable baseline
2. **Phase 2:** Add IPC push for Electron - significant UX improvement for desktop users
3. **Phase 3:** Add WebSocket/SSE for web deployment - optional, depends on infrastructure constraints

The Signal/Slot foundation makes Phase 2 and 3 natural extensions rather than rewrites.
