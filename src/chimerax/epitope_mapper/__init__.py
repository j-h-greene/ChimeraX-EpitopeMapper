from chimerax.core.commands import CmdDesc, register
from chimerax.core.toolshed import BundleAPI

class _EpitopeMapperAPI(BundleAPI):
    # Required for newer API versions
    api_version = 1 

    @classmethod
    def register_command(cls, bi, ci, logger):
        # 1. Safely import your command function from cmd.py
        from .cmd import epitope_map
        
        # 2. Define what arguments your command expects 
        desc = CmdDesc(
            synopsis='Map epitopes on protein structures.',
        )
        
        # 3. Register the command string with ChimeraX
        # Note the correct ChimeraX order: register(name, function, description)
        register(ci.name, epitope_map, desc)

    @staticmethod
    def start_tool(session, bi, ti, **kw):
        pass

    @staticmethod
    def initialize(session, bi):
        pass

    @staticmethod
    def finish(session, bi):
        pass

    @staticmethod
    def get_class(class_name):
        return None

# Export exactly ONE bundle_api instance for ChimeraX to look for
bundle_api = _EpitopeMapperAPI()