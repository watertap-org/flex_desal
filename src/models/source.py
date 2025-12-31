from pyomo.environ import Suffix

from idaes.core import declare_process_block_class
from idaes.models.unit_models.feed import FeedData

from watertap.core import InitializationMixin

from costing.source import cost_source

__author__ = "Kurban Sitterley"


@declare_process_block_class("Source")
class SourceData(InitializationMixin, FeedData):

    CONFIG = FeedData.CONFIG()

    def build(self):
        super().build()

        self.scaling_factor = Suffix(direction=Suffix.EXPORT)

    @property
    def default_costing_method(self):
        return cost_source
