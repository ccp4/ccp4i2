
import pytest
import asyncio
import time
from core.base_object.event_system import Application



@pytest.mark.asyncio
async def test_application_async_task(capsys):
    app = Application("TestApp")
    await app.start()

    async def async_task(duration: float, name: str):
        print(f"Starting async task {name}")
        for i in range(5):
            await asyncio.sleep(duration / 5)
            print(f"Task {name} progress: {(i + 1) * 20}%")
        print(f"Async task {name} completed")
        return f"Result from {name}"

    task_id = app.schedule_task(async_task, 0.5, "AsyncTask1")
    await asyncio.sleep(0.7)

    infos = app.get_task_info()
    assert any(info.name == task_id and info.state.name == "COMPLETED" for info in infos)

    captured = capsys.readouterr()
    assert "Starting async task AsyncTask1" in captured.out
    assert "Task AsyncTask1 progress: 20%" in captured.out
    assert "Async task AsyncTask1 completed" in captured.out

    await app.shutdown()


@pytest.mark.asyncio
async def test_application_sync_task(capsys):
    app = Application("TestApp")
    await app.start()

    def sync_task(duration: float, name: str):
        print(f"Starting sync task {name}")
        time.sleep(duration)
        print(f"Sync task {name} completed")
        return f"Result from {name}"

    task_id = app.schedule_task(sync_task, 0.2, "SyncTask1")
    await asyncio.sleep(0.3)

    infos = app.get_task_info()
    assert any(info.name == task_id and info.state.name == "COMPLETED" for info in infos)

    captured = capsys.readouterr()
    assert "Starting sync task SyncTask1" in captured.out
    assert "Sync task SyncTask1 completed" in captured.out

    await app.shutdown()
