import numpy as np
from numpy import exp
from numpy.linalg import inv

class EP(object):
    # constants
    BATCH_SIZE_NO = None;
    eps = np.finfo(float).eps;
    
    # variables
    constellation = None;
    constellation_len = 0;
    es = 1;                                 # constellation average power
    beta = 0.2;                             # the percentage for taking values from the previous iteration 
    epsilon = 2*5e-7;                       # the default minimal variance is 2*5e-7 (5e-7 for real only)
    l = 10;                                 # maximal iteration
    early_stop = False;
    early_stop_min_diff = 2e-4;             # the elemental minimal difference for mean and variance in Bayesian estimation (1e-4 for real only)
    batch_size = BATCH_SIZE_NO;
    
    '''
    init
    @constellation:             the constellation, a vector
    @beta:                      the percentage for taking values from the previous iteration
    @epsilon:                   the default minimal variance
    @l:                         maximal iteration
    @early_stop:                whether stop early
    @early_stop_min_diff:       the elemental minimal difference in mean and variance
    @batch_size:                the batch size
    '''
    def __init__(self, constellation, *, beta=None, epsilon=None, l=None, early_stop=None, early_stop_min_diff=None, batch_size=None):
        # inputs
        constellation = np.asarray(constellation).squeeze();
        if constellation.ndim != 1:
            raise Exception("The constellation must be a vector.");
        else:
            self.constellation = constellation;
            self.constellation_len = len(constellation);
        # inputs - constellation average power
        self.es = np.sum(self.constellation**2)/self.constellation_len;
        # optional iputs
        if beta is not None:
            self.beta = np.asarray(beta).squeeze();
        if epsilon is not None:
            self.epsilon = np.asarray(epsilon).squeeze();
        if l is not None:
            self.l = np.asarray(l).squeeze();
        if early_stop is not None:
            self.early_stop = early_stop;
        if early_stop_min_diff is not None:
            self.early_stop_min_diff = np.asarray(early_stop_min_diff).squeeze();
        if batch_size is not None:
            self.batch_size = batch_size;
        
    '''
    detect
    @y:           the received signal
    @H:           the channel matrix
    @No:          the noise (linear) power
    @sym_map:     false by default. If true, the output will be mapped to the constellation
    '''
    def detect(self, y, H, No, *, sym_map = False):
        # input check
        # input check - to numpy
        y = self.squeeze(np.asarray(y));
        H = np.asarray(H);
        No = np.squeeze(np.asarray(No));
        # input check - batch_size
        if self.batch_size != self.BATCH_SIZE_NO:
            if No.ndim == 0:
                No = np.tile(No, self.batch_size);
            if self.batch_size != y.shape[0] or self.batch_size != H.shape[0] or self.batch_size != No.shape[0]:
                raise Exception("The batch size (1st dimension) is not uniform for y, H and No.");
        # input check - dimension
        if not self.isvector(y):
            raise Exception("The received signal must be a vector.");
        if not self.ismatrix(H):
            raise Exception("The channel must be a matrix.");
        y_num = H.shape[-2];
        x_num = H.shape[-1];
        if y_num != y.shape[-1]:
            raise Exception("The channel row number does not equal to the signal number.");
        if y_num < x_num:
            raise Exception("The channel is a correlated channel.");
        if self.batch_size == self.BATCH_SIZE_NO and No.ndim > 0 or self.batch_size != self.BATCH_SIZE_NO and No.ndim > 1: 
            raise Exception("The noise power must be a scalar or a vector (only when using batch).");
        
        # reshape the input into correct shape
        # y
        y = np.expand_dims(y, -1);
        # No
        if self.batch_size != self.BATCH_SIZE_NO:
            No = np.expand_dims(No, (-1,-2));
        
        # constant values
        Ht = np.moveaxis(H, -1, -2).conj();
        Hty = Ht @ y;
        HtH = Ht @ H;
        
        # iterative detection
        gamma = self.zeros(x_num, 1);
        lamda = (1/self.es)*self.ones(x_num, 1);
        if self.early_stop:
            lmmse_mu_prev = self.zeros(x_num, 1);
            lmmse_Sigma_prev = self.zeros(x_num, x_num);
        for iter_id in range(self.l):
                # linear MMSE
                lmmse_Sigma  = inv(HtH + No*self.diag(lamda));  # covariance matrix
                lmmse_mu = lmmse_Sigma @ (Hty + No*gamma);        # mean vector
                
                # early stop
                if self.early_stop:
                    if iter_id > 0:
                        if self.batch_size == self.BATCH_SIZE_NO:
                            if np.sum(abs(lmmse_mu_prev - lmmse_mu)< self.early_stop_min_diff, axis=None) == x_num and np.sum(abs(lmmse_Sigma_prev-lmmse_Sigma)<self.early_stop_min_diff, axis=None)==x_num*x_num:
                                break;
                        else:
                            if np.sum(abs(lmmse_mu_prev - lmmse_mu)< self.early_stop_min_diff, axis=None) == self.batch_size*x_num and np.sum(abs(lmmse_Sigma_prev-lmmse_Sigma)<self.early_stop_min_diff, axis=None) == self.batch_size*x_num*x_num:
                                break;
                    # early stop - store mean and covariance
                    lmmse_mu_prev = lmmse_mu;
                    lmmse_Sigma_prev = lmmse_Sigma;
                
                # cavity marginal
                sigmai2 = No*self.diag(np.real(lmmse_Sigma));
                hi2 = sigmai2/(1-sigmai2*lamda);
                hi2 = self.max(hi2, self.eps);
                ti = hi2*(lmmse_mu/sigmai2 - gamma);
                
                # BSE - Estimate P(x|y) using Gaussian distribution
                pxyPdfExpPower = -1/(2*hi2)*abs(self.repmat(ti, 1, self.constellation_len) - self.repmat(self.constellation, x_num, 1))**2;
                # BSE - make every row the max power is 0
                #     - max only consider the real part
                pxypdfExpNormPower = pxyPdfExpPower - np.expand_dims(pxyPdfExpPower.max(axis=-1), axis=-1);
                #pxypdfExpNormPower = pxyPdfExpPower - np.expand_dims(self.max(pxyPdfExpPower, axis=-1), axis=-1);
                pxyPdf = exp(pxypdfExpNormPower);
                # BSE - Calculate the coefficient of every possible x to make the sum of all
                pxyPdfCoeff = np.expand_dims(1./np.sum(pxyPdf, axis=-1), -1);
                pxyPdfCoeff = self.repmat(pxyPdfCoeff, 1, self.constellation_len);
                # BSE - PDF normalisation
                pxyPdfNorm = pxyPdfCoeff*pxyPdf;
                # BSE - calculate the mean and variance
                mu_pi = np.expand_dims(np.sum(pxyPdfNorm*self.constellation, axis=-1), -1);
                mu_pi_mat = self.repmat(mu_pi, 1, self.constellation_len);
                sigma_pi = np.expand_dims(np.sum(abs(mu_pi_mat - self.constellation)**2*pxyPdfNorm, axis=-1), -1);
                sigma_pi = self.max(sigma_pi, self.epsilon);
                
                # update
                lamda_new = (hi2-sigma_pi)/sigma_pi/hi2;
                gamma_new = (mu_pi*hi2-sigma_pi*ti)/sigma_pi/hi2;
                # update - avoid Negative Index
                lamda_new_neg_idx = lamda_new < 0;
                lamda_new[lamda_new_neg_idx] = lamda[lamda_new_neg_idx];
                gamma_new[lamda_new_neg_idx] = gamma[lamda_new_neg_idx];
                # update - convex combination
                lamda = (1-self.beta)*lamda + self.beta*lamda_new;
                gamma = (1-self.beta)*gamma + self.beta*gamma_new;
        x = mu_pi.squeeze(-1);
        if sym_map:
            x = self.symmap(x);
        return x;
                
        
    '''
    symbol mapping (hard)
    '''
    def symmap(self, syms):
        syms = np.asarray(syms);
        if not self.isvector(syms):
            raise Exception("Symbols must be into a vector form to map.");
        syms_len = syms.shape[-1];
        syms_mat = self.repmat(np.expand_dims(syms, -1), 1, self.constellation_len);
        constel_mat = self.repmat(self.constellation, syms_len, 1);
        syms_dis = abs(syms_mat - constel_mat)**2;
        syms_dis_min_idx = syms_dis.argmin(axis=-1);
        
        return np.take(self.constellation, syms_dis_min_idx);
    
    
    ##########################################################################
    # Functions uniform with non-batch and batch
    ##########################################################################
    '''
    check input is a vector like [(batch_size), n],  [(batch_size), n, 1] or [(batch_size), 1, n] 
    '''
    def isvector(self, mat):
        mat = np.asarray(mat);
        if self.batch_size == self.BATCH_SIZE_NO:
            return mat.ndim == 1 or mat.ndim == 2 and (mat.shape[-2] == 1 or mat.shape[-1] == 1);
        else:
            if mat.shape[0] != self.batch_size:
                raise Exception("The input does not has the required batch size.");
            else:
                return mat.ndim == 2 or mat.ndim == 3 and (mat.shape[-2] == 1 or mat.shape[-1] == 1);
    
    '''
    check input is a matrix like [(batch_size), n. m]
    '''
    def ismatrix(self, mat):
        mat = np.asarray(mat);
        if self.batch_size == self.BATCH_SIZE_NO:
            return mat.ndim == 2 and mat.shape[-2] > 1 and mat.shape[-1] > 1;
        else:
            if mat.shape[0] != self.batch_size:
                raise Exception("The input does not has the required batch size.");
            else:
                return mat.ndim == 3 and mat.shape[-2] > 1 and mat.shape[-1] > 1;
    
    
    '''
    return the maximum value of a matrix or the maximu value of two matrices (for complex value, we compare the magnitude)
    @mat: the matrix to 
    '''
    def max(self, mat, *args, axis=-1):
        out = None;
        mat = np.asarray(mat);
        mat2 = None;
        if len(args) > 0:
            if isinstance(args[0], float) or isinstance(args[0], int):
                mat2 = args[0];
            elif args[0].ndim == 0 or mat.shape == args[0].shape:
                mat2 = args[0].astype(mat.dtype);
            else:
                raise Exception("The input two matrices must have the same shape.");
        # non-complex value
        if mat.dtype != np.csingle and mat.dtype != np.complex_:
            if len(args) == 0:
                out = mat.max(axis=axis);
            else:
                out = np.where(mat>mat2, mat, mat2);
        # complex value
        else:
            if len(args) == 0:
                out = np.take_along_axis(mat, abs(mat).argmax(axis=axis, keepdims=True), axis).squeeze(axis);
            else:
                out = np.where(abs(mat)>abs(mat2), mat, mat2);
        return out;
    
    '''
    squeeze redundant dimension except the batch size
    '''
    def squeeze(self, mat):
        out = mat.copy();
        out = np.squeeze(out);
        if self.batch_size == 1:
            out = np.expand_dims(out, 0);
        return out;
    
    '''
    return an identity matrix
    @size: a tuple of the shape
    '''
    def eye(self, size):
        out = np.eye(size);
        if self.batch_size is not self.BATCH_SIZE_NO:
            out = np.tile(out,(self.batch_size, 1, 1));
        return out;
    
    '''
    generate a matrix of all zeros
    @order: 'C': this function only create given dimensions; 'F': create the dimensions as matlab (2D at least)
    '''
    def zeros(self, nrow, *args, order='C'):
        out = None;
        if order == 'F':
            ncol = nrow;
            if len(args) >= 1:
                ncol = args[0];
            out = np.zeros((nrow, ncol)) if self.batch_size == self.BATCH_SIZE_NO else np.zeros((self.batch_size, nrow, ncol));
        elif order == 'C':
            zeros_shape = list(args);
            zeros_shape.insert(0, nrow);
            if self.batch_size != self.BATCH_SIZE_NO:
                zeros_shape.insert(0, self.batch_size);
            out = np.zeros(zeros_shape);
        return out;
    
    '''
    generate a matrix full of ones
    @order: 'C': this function only create given dimensions; 'F': create the dimensions as matlab (2D at least)
    '''
    def ones(self, nrow, *args, order='C'):
        shape = [];
        # format the shape for dimensions
        if order == 'F':
            ncol = nrow;
            if len(args) >= 1:
                ncol = args[0];
            shape.append(nrow);
            shape.append(ncol);
        elif order == 'C':
            shape = list(args);
            shape.insert(0, nrow);
        else:
            raise Exception("The order is illegal.");
        # format the shape for batch_size
        if self.batch_size != self.BATCH_SIZE_NO:
            shape.insert(0, self.batch_size);
        return np.ones(shape);
        
    
    '''
    generate a matrix based on its diag or get a diagonal matrix from its vector
    @mat:   a vector as [(batch_size), n, 1] or [(batch_size), 1, n]; if n == 1, 
            it will be taken as [(batch_size), n, 1].
            Or a square matrix [(batch_size), n, n]
            
    '''
    def diag(self, mat):
        # remove redundant dimensions
        if self.batch_size == self.BATCH_SIZE_NO:
            if mat.ndim >= 2:
                if mat.shape[-1] == 1:
                    mat = mat.squeeze(-1);
                elif mat.shape[-2] == 1:
                    mat = mat.squeeze(-2);
        else:
            if mat.ndim >= 3:
                if mat.shape[-1] == 1:
                    mat = mat.squeeze(-1);
                elif mat.shape[-2] == 1:
                    mat = mat.squeeze(-2);
        # generate output
        out = None;
        # diag_vec_len = diag_vec.shape[-1];
        if self.batch_size is self.BATCH_SIZE_NO:
            out = np.diag(mat);
        else:
            # np.zeros only take real numbers by default, here we need to put complex value into it
            out = [];
            # create output
            for batch_id in range(self.batch_size):
                out.append(np.diag(mat[batch_id, ...]));
            out = np.asarray(out);
        # add extra dimension in the end if we want a vector
        if self.ismatrix(mat):
            out = np.expand_dims(out, axis=-1);
        return out;
    
    '''
    repeat the matrix in the given dimension (as matlab)
    '''
    def repmat(self, mat, nrow, *args):
        out = None;
        ncol = args[0] if len(args) >= 1 else nrow;
        out = np.tile(mat, (nrow, ncol)) if self.batch_size == self.BATCH_SIZE_NO else np.tile(mat, (1, nrow, ncol));
        return out;