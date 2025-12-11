import giant.logs as lg
from giant.mulch.statistics import ExtractAndClassifyDatasetStatistics
from giant.mulch.statistics.scaling import (
    ClassifyScalingStatistics,
    ExtractScalingStatistics,
)
from giant.mulch.statistics.xray import (
    ClassifyWilsonStatistics,
    ExtractBasicXrayStatistics,
    ExtractWilsonStatistics,
)

logger = lg.getLogger(__name__)


class PanddaExtractDatasetStatistics(ExtractAndClassifyDatasetStatistics):

    def __init__(
        self,
        max_scaling_z_score,
        max_wilson_z_score,
    ):

        extracters = [
            ExtractBasicXrayStatistics(),
            ExtractWilsonStatistics(),
            ExtractScalingStatistics(),
        ]

        classifiers = [
            ClassifyWilsonStatistics(
                outlier_partition="not_train",
                max_z_score=max_wilson_z_score,
            ),
            ClassifyScalingStatistics(
                outlier_partition="not_train",
                max_z_score=max_scaling_z_score,
            ),
        ]

        super(PanddaExtractDatasetStatistics, self).__init__(
            extracters=extracters,
            classifiers=classifiers,
        )
