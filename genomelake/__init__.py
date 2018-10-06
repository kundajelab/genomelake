from __future__ import absolute_import

__version__ = '0.1.6'

# Setup logging
import logging

log_formatter = \
    logging.Formatter('%(levelname)s:%(asctime)s:%(name)s] %(message)s')
_logger = logging.getLogger('genomelake')
_handler = logging.StreamHandler()
_handler.setLevel(logging.DEBUG)
_handler.setFormatter(log_formatter)
_logger.setLevel(logging.DEBUG)
_logger.addHandler(_handler)
