"""Train-only protein/context representations with explicit solvent interactions."""
import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


def products(left,right):
    return (left[:,:,None]*right[:,None,:]).reshape(len(left),-1)


class InteractionFeatures(TransformerMixin,BaseEstimator):
    def __init__(self,conditions,protein_columns,pca_components=None):
        self.conditions=conditions
        self.protein_columns=protein_columns
        self.pca_components=pca_components

    def fit(self,frame,y=None):
        self.context_=Pipeline([
            ('impute',SimpleImputer(strategy='median',add_indicator=True,keep_empty_features=True)),
            ('scale',StandardScaler())]).fit(frame[self.conditions])
        self.solvents_=OneHotEncoder(handle_unknown='ignore',sparse_output=False).fit(frame[['solvent_name']])
        if self.protein_columns:
            unique=frame.drop_duplicates('sequence_id')
            values=unique[self.protein_columns].to_numpy()
            if not np.isfinite(values).all():
                raise ValueError('Protein descriptors must be finite')
            self.pca_=None
            if self.pca_components is not None:
                self.pca_=PCA(n_components=min(self.pca_components,len(unique)-1,len(self.protein_columns)),
                              svd_solver='randomized',random_state=11)
                values=self.pca_.fit_transform(values)
            self.protein_scaler_=StandardScaler().fit(values)
        return self

    def transform(self,frame):
        context=self.context_.transform(frame[self.conditions])
        solvents=self.solvents_.transform(frame[['solvent_name']])
        pieces=[context,solvents,products(solvents,context)]
        if self.protein_columns:
            values=frame[self.protein_columns].to_numpy()
            if self.pca_ is not None:
                values=self.pca_.transform(values)
            protein=self.protein_scaler_.transform(values)
            pieces.extend([protein,products(protein,solvents),products(protein,context)])
        return np.concatenate(pieces,axis=1)
