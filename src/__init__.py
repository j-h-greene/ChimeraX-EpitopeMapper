from chimerax.core.commands import CmdDesc, CmdDescArg, CmdDescKeyword, register
from chimerax.core.toolshed import BundleAPI

class _EpitopeMapperAPI(BundleAPI):
    
    # Required for newer API versions
    api_version = 1 

    @staticmethod
    def register_command(bi, ci, logger):
        from .cmd import epitope_map
        
        # We define the command description here (moved from your old code)
        desc = CmdDesc(
            required=[('atoms', CmdDescArg(str, help="The atoms to evaluate (e.g., :1)"))],
            keyword=[
                ('cutoff', CmdDescKeyword(float, help="Distance cutoff for spatial clustering")),
            ],
            synopsis='Map epitopes on protein structures.',
        )
        register(ci.name, desc, epitope_map, logger=logger)

    @staticmethod
    def start_tool(session, bi, ti, **kw):
        pass

    # --- NEW LIFECYCLE METHODS ---

    @staticmethod
    def initialize(session, bi):
        # Called each time ChimeraX starts if custom initialization is set.
        pass

    @staticmethod
    def finish(session, bi):
        # Called when the bundle is deinitialized/unloaded. 
        # Good for cleaning up global variables or temporary files.
        pass

    @staticmethod
    def get_class(class_name):
        # Called by session code to get a class from the bundle that was saved in a session.
        # If you later create a UI tool (e.g., EpitopeMapperUI), return it here.
        return None

bundle_api = _EpitopeMapperAPI()