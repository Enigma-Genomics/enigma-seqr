import luigi

from loading_pipeline.lib.misc.callset_identity import guard_callset
from loading_pipeline.lib.tasks.base.base_loading_run_params import (
    BaseLoadingRunParams,
)
from loading_pipeline.lib.tasks.dataproc.base_run_job_on_dataproc import (
    BaseRunJobOnDataprocTask,
)
from loading_pipeline.lib.tasks.run_pipeline import RunPipelineTask
from loading_pipeline.lib.tasks.write_existing_variants_parquet import (
    WriteExistingVariantsParquetTask,
)


@guard_callset
@luigi.util.inherits(BaseLoadingRunParams)
class RunPipelineOnDataprocTask(BaseRunJobOnDataprocTask):
    attempt_id = luigi.IntParameter()

    @property
    def task(self) -> luigi.Task:
        return RunPipelineTask

    def requires(self) -> [luigi.Task]:
        return [self.clone(WriteExistingVariantsParquetTask), *super().requires()]
