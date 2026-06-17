from chimerax.core.commands import CmdDesc, CmdDescArg, CmdDescKeyword, register
from chimerax.core.toolshed import BundleAPI

class _EpitopeMapperAPI(BundleAPI):
    @staticmethod
    def register_command(command_name, logger):
        from .epitope_mapper import EpitopeMapper
        desc = CmdDesc(
            required=[('model', CmdDescArg(str, help="The model ID"))],
            keyword=[
                ('output', CmdDescKeyword(str, help="Output file path")),
            ],
            synopsis='Map epitopes on protein structures.',
        )
        register(command_name, desc, EpitopeMapper, logger=logger)

    @staticmethod
    def start_tool(session, tool_name, **kw):
        pass

bundle_api = _EpitopeMapperAPI()