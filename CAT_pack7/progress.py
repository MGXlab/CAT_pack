from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

from rich import filesize
from rich.console import Console
from rich.progress import BarColumn, Progress, ProgressColumn, Task, TextColumn, TimeElapsedColumn
from rich.text import Text

from .pipeline import Step
from .utils.logging import Status


class DetailsColumn(ProgressColumn):
    def render(self, task: Task) -> Text:
        state, unit = task.fields["state"], task.fields["step"].unit
        parts = [str(state)]
        if unit == "bytes":
            total = filesize.decimal(task.total) if task.total is not None else "?"
            parts.append(f"{filesize.decimal(task.completed)} / {total}")
        if unit and task.speed is not None and state is Status.RUNNING:
            speed = filesize.decimal(task.speed) if unit == "bytes" else f"{task.speed:,.1f} {unit}"
            parts.append(f"{speed}/s")
        return Text(" | ".join(parts))


class GroupedProgress(Progress):
    def get(self):
        rows = []
        for task in self.tasks:
            step = task.fields["step"]
            if step.parent is None:
                rows.append(task)
                rows.extend(subtask for subtask in self.tasks if subtask.fields["step"].parent == step.name)
        yield self.make_tasks_table(rows)


@dataclass
class ProgressTracker:
    steps: list[Step]
    console: Console
    quiet: bool = False
    progress: Progress = field(init=False)
    tasks: dict[str, Task] = field(default_factory=dict, init=False)
    active: Task | None = field(default=None, init=False)
    last_update: float = field(default=0.0, init=False)

    def __post_init__(self)-> None:
        self.progress = GroupedProgress(
            TextColumn("{task.description}", markup=False),
            BarColumn(bar_width=20),
            DetailsColumn(),
            TimeElapsedColumn(),
            console=self.console, disable=self.quiet,
            refresh_per_second=1, speed_estimate_period=60,
        )

        for step in self.steps:
            if step.parent and step.parent not in self.tasks:
                self.add(Step(step.parent))
            self.add(step, visible=step.parent is None)

        for name, task in self.tasks.items():
            if subprocesses := self.subprocesses(name):
                self.progress.update(task.id, total=len(subprocesses))

    def add(self, step: Step, total=1, visible=True) -> None:
        label = step.display_name or step.name
        self.progress.add_task(f"  {label}" if step.parent else label, total=total,
                               start=False, visible=visible, step=step, state=Status.WAITING)
        self.tasks[step.name] = self.progress.tasks[-1]

    def subprocesses(self, name: str) -> list[Task]:
        return [task for task in self.tasks.values() if task.fields["step"].parent == name]

    def set(self, task: Task, status: Status, completed=None, total=None):
        if status is Status.RUNNING: self.progress.start_task(task.id)
        elif status is not Status.WAITING: self.progress.stop_task(task.id)
        self.progress.update(task.id, state=status, completed=completed, total=total)

    def update(self, name: str, status: Status, completed=0, total=None):
        task = self.tasks[name]
        changed = status is not task.fields["state"]
        done = (Status.COMPLETE, Status.SUPPLIED, Status.SKIPPED)
        count = completed if total is not None else None

        if status is not Status.RUNNING:
            count, total = task.total if status in done else None, None
        self.set(task, status, count, total)

        if parent := task.fields["step"].parent:
            subprocesses = self.subprocesses(parent)
            if not task.visible and status not in (Status.WAITING, Status.SKIPPED):
                for subprocess in subprocesses:
                    self.progress.update(subprocess.id, visible=True)
            proc_states = [subprocess.fields["state"] for subprocess in subprocesses]
            count = sum(state in done for state in proc_states)
            count += sum(subprocess.percentage / 100 for subprocess in subprocesses
                         if subprocess.fields["state"] is Status.RUNNING)
            if Status.FAILED in proc_states or Status.CANCELLED in proc_states:
                status = Status.FAILED if Status.FAILED in proc_states else Status.CANCELLED
            elif all(state in done for state in proc_states):
                status = Status.SKIPPED if all(state is Status.SKIPPED for state in proc_states) else Status.COMPLETE
            else:
                status = Status.WAITING if all(state is Status.WAITING for state in proc_states) else Status.RUNNING
            self.set(self.tasks[parent], status, count)

        if changed: self.progress.refresh()

    def queue(self, paths: list[Path]):
        for path in paths:
            self.add(Step(str(path), parent="Downloading", display_name=path.name, unit="bytes"), total=None)

    def report_download(self, path: Path, completed: int, total: int | None):
        task = self.tasks[str(path)]
        now = self.progress.get_time()
        if task is not self.active:
            self.finish(Status.COMPLETE)
            self.active, self.last_update = task, now - 0.5
        if completed > 0 and (total is None or completed < total) and now - self.last_update < 0.5:
            return
        self.set(task, Status.RUNNING, completed, total)
        self.last_update = now

    def finish(self, status: Status):
        if self.active is not None:
            self.set(self.active, Status.COMPLETE if self.active.finished else status)

    @contextmanager
    def display(self):
        with self.progress:
            status = Status.COMPLETE
            try:
                yield
            except KeyboardInterrupt:
                status = Status.CANCELLED
                raise
            except Exception:
                status = Status.FAILED
                raise
            finally:
                self.finish(status)
