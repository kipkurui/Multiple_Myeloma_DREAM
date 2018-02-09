import synapseclient
import os
from synapseclient import Activity
from synapseclient import Entity, Project, Folder, File, Link
from synapseclient import Evaluation, Submission, SubmissionStatus
from synapseclient import Wiki

syn = synapseclient.Synapse()
syn.login('kurui','AchameJeso@17')

dest_dir = '/mnt/lustre/users/ckibet/MM_DREAM/Data/'

import synapseutils
files = synapseutils.syncFromSynapse(syn, 'syn7222203',path=dest_dir)