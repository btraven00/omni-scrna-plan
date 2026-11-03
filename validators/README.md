# validators

By convention, we drop ad-hoc validators that are executed for every file produced by a stage.

* Every validation script receives a single file path.
* Validators are dropped at `STAGE/OUTPUT_ID` (stage and output ids from the plan).
* Use as few dependencies as possible. Benchmarker adds them to the `validation` pixi environment in `pixi.toml'
