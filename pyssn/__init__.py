
__all__ = []
__version__ = '0.2.2'

from utils.Config import _Config
config = _Config()
log_ = config.log_


log_.message('Starting pySSN.', calling = 'pySSN init')

from core.spectrum import read_data, spectrum
from utils.physics import CST
