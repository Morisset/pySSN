#!/usr/bin/env python
if __name__ == '__main__':
    import pyssn
    import sys

    if pyssn.config.INSTALLED['Qt4']:
        from pyssn.qt.pyssn_qt import main
        sys.exit(main())
    else:
        from pyssn.core.spectrum import main
        main()
