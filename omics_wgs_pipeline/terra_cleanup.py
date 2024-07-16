# %%

from omics_wgs_pipeline.terra import call_firecloud_api
import pandas as pd
import requests
import re
from firecloud import api as firecloud_api
from google.cloud import storage
import google.auth.transport.requests
import google.oauth2.id_token
import os

os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/Users/golden/Documents/secrets/omics/terra_omics_google_auth.json"


abort_suff = '/api/workflows/{version}/{id}/abort'
fcrest_pref = 'https://api.firecloud.org'

def rest_api_abort(wfid):
    token = get_gcp_oidc_token()
    api_url = fcrest_pref + abort_suff.format(version='v1',id=wfid)
    command = {}
    response = requests.post(api_url,headers={"Authorization": f"Bearer {token}"})
    return response


def get_gcp_oidc_token() -> str:
    """
    Get a GCP OIDC token (ID token) for current credentials.

    :return: the auth/bearer token
    """

    auth_req = google.auth.transport.requests.Request()  # type: ignore

    token = google.oauth2.id_token.fetch_id_token(  # type: ignore
        auth_req, "https://cloudfunctions.googleapis.com"
    )

    if token is None:
        raise ValueError("GCP auth token cannot be None")

    return token


# %%
client = storage.Client('depmap-omics')

bucket = client.get_bucket('fc-secure-0a4e826a-6c7e-4481-89be-42dfb76f862d')

blob_pref = 'submissions/{submission}/annotateVariants/{workflowId}/call-open_cravat/'
blob_format = 'submissions/{submission}/annotateVariants/{workflowId}/call-open_cravat/stderr'
submission_id = '1c8a752b-e599-4a25-bede-5ee97edd192f'

df_submission = call_firecloud_api(
    firecloud_api.get_submission,
    namespace='broad-firecloud-ccle',
    workspace='tcga_mutation_testing',
    submission_id=submission_id,
)

df_workflows = pd.DataFrame.from_dict(data=df_submission['workflows'])
df_workflows = df_workflows[df_workflows.status == 'Running']

match_strs = []
match_strs.append('ERROR: Could not install packages due to an OSError: [Errno 2]')
match_strs.append('ERROR: Exception')


def get_err_text(wf_id):
    blobs = list(bucket.list_blobs(prefix=blob_pref.format(submission=submission_id, workflowId=wf_id)))
    stderrs = [b for b in blobs if 'stderr' in b.name]
    if len(stderrs) > 1:
        stderrs = sorted(stderrs, key=lambda x: re.search(r'call-open_cravat\/(.+)', x.name).group(1))
        # the last element will be the first attempt, the next will be the most recent
        text = stderrs[-2].download_as_text()
    elif len(stderrs) == 0:
        return None
    else:
        text = stderrs[0].download_as_text()

    return text


stuck_wfs = []

for i, wf_id in enumerate(df_workflows.workflowId):
    try:
        err_text = get_err_text(wf_id)
    except AttributeError:
        err_text = None

    if err_text is not None and any([s in err_text for s in match_strs]):
        stuck_wfs.append(wf_id)

    if len(df_workflows) > 9 and i % (len(df_workflows) // 10) == 0:
        print(100 * i // len(df_workflows), '% done')

print(str(len(stuck_wfs))+' stalled workflows')
pd.Series(stuck_wfs).to_csv('stuck_wfs.csv',index=False)
print('aborting '+str(len(stuck_wfs))+' stalled workflows')

for i, wf in enumerate(stuck_wfs):
    rest_api_abort(wf)
    if len(stuck_wfs) > 9 and i % (len(stuck_wfs) // 10) == 0:
        print(100 * i // len(stuck_wfs), '% done')