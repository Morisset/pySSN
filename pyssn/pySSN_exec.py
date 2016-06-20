#!/usr/bin/env python

def main():
    import pyssn
    import sys
    from pyssn.utils.misc import get_parser

    parser = get_parser()
    args = parser.parse_args()
    if args.noQt4:
        pyssn.config.INSTALLED['Qt4'] = False

    if pyssn.config.INSTALLED['Qt4']:
        from pyssn.qt.pyssn_qt import main
        sys.exit(main())
    else:
        from pyssn.core.spectrum import main
        sys.exit(main())
    
if __name__ == '__main__':
    main()
    
