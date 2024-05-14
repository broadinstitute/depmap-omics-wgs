from typing import TypedDict

import pandera.typing

TypedDataFrame = pandera.typing.DataFrame


class PersistedWdl(TypedDict):
    wdl: str
    public_url: str
