import sys

from mock import Mock
from mock import patch

# uncomment if cctbx not available

# mock cctbx modules
# add them if they are not installed
# sys.modules["iotbx"] = Mock()
# sys.modules["iotbx.pdb"] = Mock()
#
# with patch('iotbx.pdb', Mock(return_value="")):
#     pass
