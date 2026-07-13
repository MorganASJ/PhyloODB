"""Analysis task implementations."""
from .busco import DownloadBuscoLibraryTask, BatchBuscoTask, BuscoTask
from .orthofinder import OrthoFinderTask
from .add_library import AddLibraryTask
from .import_custom_library import ImportCustomLibraryTask
from .blastdb import CreateProteomeBlastDB, ConstructBuscoBlastDB
from .paralog_removal import ParalogRemovalTask
from .decontamination import Decontamination
from .internal_decontamination import InternalDecontaminationTask
from .external_decontamination_check import ExternalDecontaminationCheckTask
from .external_decontamination_apply import ExternalDecontaminationApplyTask
from .trees import MafftTask, IQTreeTask, BuildBuscoTreesTask, AnnotateOrthogroupTreeTask

__all__ = [
    "DownloadBuscoLibraryTask",
    "BatchBuscoTask",
    "BuscoTask",
    "OrthoFinderTask",
    "AddLibraryTask",
    "ImportCustomLibraryTask",
    "CreateProteomeBlastDB",
    "ConstructBuscoBlastDB",
    "ParalogRemovalTask",
    "Decontamination",
    "InternalDecontaminationTask",
    "ExternalDecontaminationCheckTask",
    "ExternalDecontaminationApplyTask",
    "MafftTask",
    "IQTreeTask",
    "BuildBuscoTreesTask",
    "AnnotateOrthogroupTreeTask",
]
