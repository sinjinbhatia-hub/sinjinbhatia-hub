# Build workflow runs (for resuming after interruptions)

Script: `curriculum/tools/build_modules_workflow.js`. Resume a run with
`Workflow({scriptPath: <script>, resumeFromRunId: <id>, args: <same args>})`; completed agents replay from cache.

| Run ID | Modules |
|---|---|
| wf_6c57abdb-8b7 | verify+fix: LA-1..LA-12, LA-14, F1 |
| wf_08ca99dc-1c9 | write: F4, F5; verify: F2, F3, F6 |
| wf_fda2e292-c61 | write: LA-13, CM-1, CM-2 |
| wf_476c56eb-6ca | write: CM-3, CM-4, CM-5 |
| wf_5f4c5a5d-78d | write: CM-6, CM-7, CM-8 |
| wf_6814d784-9c6 | write: EM-1..EM-4 |
| wf_9264db32-034 | write: QM-1..QM-4 |
| wf_9f925196-492 | write: SM-1..SM-4 |
