__all__ = []
from .version import __version__

from utils.Config import _Config
config = _Config()
log_ = config.log_
log_.message('Starting pySSN.', calling = 'pySSN init')

from core.spectrum import read_data, spectrum
from utils.physics import CST
if config.INSTALLED['Qt4']:
    from qt import pyssn_qt
