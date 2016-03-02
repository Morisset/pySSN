"""
pySSN is available under the GNU licence providing you cite the developpers names:

    Ch. Morisset (Instituto de Astronomia, Universidad Nacional Autonoma de Mexico)

    D. Pequignot (Meudon Observatory, France)
"""

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
