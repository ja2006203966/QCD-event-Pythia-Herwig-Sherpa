import tensorflow as tf
from tensorflow import keras
from mygraph_function import supp

# A shape is (N, P_A, C), B shape is (N, P_B, C)
# D shape is (N, P_A, P_B)
def batch_distance_matrix_general(A, B):
    with tf.name_scope('dmat'):
        r_A = tf.reduce_sum(A * A, axis=2, keepdims=True)
        r_B = tf.reduce_sum(B * B, axis=2, keepdims=True)
        m = tf.matmul(A, tf.transpose(B, perm=(0, 2, 1)))
        D = r_A - 2 * m + tf.transpose(r_B, perm=(0, 2, 1))
        return D
    
    
    

def knn(num_points, k, topk_indices, features):
    # topk_indices: (N, P, K)
    # features: (N, P, C)
    with tf.name_scope('knn'):
        queries_shape = tf.shape(features)
        batch_size = queries_shape[0]
        batch_indices = tf.tile(tf.reshape(tf.range(batch_size), (-1, 1, 1, 1)), (1, num_points, k, 1))
        indices = tf.concat([batch_indices, tf.expand_dims(topk_indices, axis=3)], axis=3)  # (N, P, K, 2)
        return tf.gather_nd(features, indices)

##====================================================================================================

## write the gpu function
def myweight(A,B,w=-2, actype = 'max'):
    with tf.name_scope('myweight'):    
        b = tf.tile(A,[1,1,100])
        b = tf.expand_dims(b,-1)
        c = tf.transpose(B, perm= (0, 2, 1))
        c = tf.tile(c,[1,100,1])
        c = tf.expand_dims(c,-1)
        W = tf.concat([b,c] ,-1)
        if actype == 'max':
            W = tf.reduce_max(W, axis=-1)
        if actype == 'min':
            W = tf.reduce_min(W, axis=-1)
        if actype == 'mean':
            W = tf.reduce_mean(W, axis=-1)
        if w < 0 :
            W = 1/W
            w = -w
        return W**w
        

##=======================================================================================================
def edge_conv(points, ptw, features, num_points, K, channels,  w=-1, actype = 'max',with_bn=True, activation='relu', pooling='average', name='edgeconv'):
    """EdgeConv
    Args:
        K: int, number of neighbors
        in_channels: # of input channels
        channels: tuple of output channels
        pooling: pooling method ('max' or 'average')
    Inputs:
        points: (N, P, C_p)
        features: (N, P, C_0)
    Returns:
        transformed points: (N, P, C_out), C_out = channels[-1]
    """

    with tf.name_scope('edgeconv'):

        # distance
        D = batch_distance_matrix_general(points, points)  # (N, P, P)
#         if points.shape[0] is not None:
        W = myweight(ptw,ptw,w=w*2 , actype = actype)
        D = D*W
        _, indices = tf.nn.top_k(-D, k=K + 1)  # (N, P, K+1)
        indices = indices[:, :, 1:]  # (N, P, K)
# =======================================================
#         K2 = K//3*2
#         K3 = K//3
#         _, indices2 = tf.nn.top_k(-D, k=K2 + 1)  # (N, P, K+1)
#         indices2 = indices2[:, :, 1:] 
#         _, indices3 = tf.nn.top_k(-D, k=K3 + 1)  # (N, P, K+1)
#         indices3 = indices3[:, :, 1:] 
# ====================================================
        
#         indices = tf.constant(neighborhood(points, k=K, jetpt = 1 , w=-1, actype='min')  )
        fts = features
        knn_fts = knn(num_points, K, indices, fts)  # (N, P, K, C)
        knn_fts_center = tf.tile(tf.expand_dims(fts, axis=2), (1, 1, K, 1))  # (N, P, K, C)
        knn_fts = tf.concat([knn_fts_center, tf.subtract(knn_fts, knn_fts_center)], axis=-1)  # (N, P, K, 2*C)
#==============================================================================================
#         fts2 = features
#         knn_fts2 = knn(num_points, K2, indices2, fts2)  # (N, P, K, C)
#         knn_fts_center2 = tf.tile(tf.expand_dims(fts2, axis=2), (1, 1, K2, 1))  # (N, P, K, C)
#         knn_fts2 = tf.concat([knn_fts_center2, tf.subtract(knn_fts2, knn_fts_center2)], axis=-1)  # (N, P, K, 2*C)
#         fts3 = features
#         knn_fts3 = knn(num_points, K3, indices3, fts3)  # (N, P, K, C)
#         knn_fts_center3 = tf.tile(tf.expand_dims(fts3, axis=2), (1, 1, K3, 1))  # (N, P, K, C)
#         knn_fts3 = tf.concat([knn_fts_center3, tf.subtract(knn_fts3, knn_fts_center3)], axis=-1)  # (N, P, K, 2*C)
#========================================================================================
#         knn_fts = tf.concat([knn_fts_center, knn_fts], axis=-1)  # (N, P, K, 2*C)

        x = knn_fts
#         x2 = knn_fts2
#         x3 = knn_fts3
        for idx, channel in enumerate(channels):
            x = keras.layers.Conv2D(channel, kernel_size=(1, 1), strides=1, data_format='channels_last',
                                    use_bias=False if with_bn else True, kernel_initializer='glorot_normal', name='%s_conv%d' % (name, idx))(x)
            if with_bn:
                x = keras.layers.BatchNormalization(name='%s_bn%d' % (name, idx))(x)
            if activation:
                x = keras.layers.Activation(activation, name='%s_act%d' % (name, idx))(x)
# ===========================================================================================================
#         for idx, channel in enumerate(channels):
#             x2 = keras.layers.Conv2D(channel, kernel_size=(1, 1), strides=1, data_format='channels_last',
#                                     use_bias=False if with_bn else True, kernel_initializer='glorot_normal', name='%s_2conv%d' % (name, idx))(x2)
#             if with_bn:
#                 x2 = keras.layers.BatchNormalization(name='%s_2bn%d' % (name, idx))(x2)
#             if activation:
#                 x2 = keras.layers.Activation(activation, name='%s_2act%d' % (name, idx))(x2)
#         for idx, channel in enumerate(channels):
#             x3 = keras.layers.Conv2D(channel, kernel_size=(1, 1), strides=1, data_format='channels_last',
#                                     use_bias=False if with_bn else True, kernel_initializer='glorot_normal', name='%s_3conv%d' % (name, idx))(x3)
#             if with_bn:
#                 x3 = keras.layers.BatchNormalization(name='%s_3bn%d' % (name, idx))(x3)
#             if activation:
#                 x3 = keras.layers.Activation(activation, name='%s_3act%d' % (name, idx))(x3)

# ===================================================================
#         x = tf.concat([x,x2,x3], axis=2)
#         x = tf.math.reduce_max(x, axis = -1)
        if pooling == 'max':
            fts = tf.reduce_max(x, axis=2)  # (N, P, C')
        else:
            fts = tf.reduce_mean(x, axis=2)  # (N, P, C')

        # shortcut
        sc = keras.layers.Conv2D(channels[-1], kernel_size=(1, 1), strides=1, data_format='channels_last',
                                 use_bias=False if with_bn else True, kernel_initializer='glorot_normal', name='%s_sc_conv' % name)(tf.expand_dims(features, axis=2))
        if with_bn:
            sc = keras.layers.BatchNormalization(name='%s_sc_bn' % name)(sc)
        sc = tf.squeeze(sc, axis=2)

        if activation:
            return keras.layers.Activation(activation, name='%s_sc_act' % name)(sc + fts)  # (N, P, C')
        else:
            return sc + fts


def _particle_net_base(points, features=None, mask=None, setting=None, w=-1, actype = 'max', name='particle_net'):
    # points : (N, P, C_coord)
    # features:  (N, P, C_features), optional
    # mask: (N, P, 1), optinal

    
    with tf.name_scope(name):
        if features is None:
            features = points[:,:,:2]

        if mask is not None:
            mask = tf.cast(tf.not_equal(mask, 0), dtype='float32')  # 1 if valid
            coord_shift = tf.multiply(999., tf.cast(tf.equal(mask, 0), dtype='float32'))  # make non-valid positions to 99

        fts = tf.squeeze(keras.layers.BatchNormalization(name='%s_fts_bn' % name)(tf.expand_dims(features, axis=2)), axis=2)
        
        
        for layer_idx, layer_param in enumerate(setting.conv_params):
            K, channels = layer_param
            pts =  points[:, :, :2]#tf.add(coord_shift, points[:, :, :2]) # if layer_idx == 0 else tf.add(coord_shift, fts)
#             mcgcnn = MCECBlock(pts, tf.expand_dims(fts[:, :, -1], axis = -1), fts, setting.num_points, K, channels, with_bn=True, w=w, actype = actype, activation='relu',pooling=setting.conv_pooling, name='%s_%s%d' % (name, 'EdgeConv', layer_idx))
#             fts = mcgcnn.call()
            fts = edge_conv(pts, tf.expand_dims(fts[:, :, -1], axis = -1), fts, setting.num_points, K, channels, with_bn=True, w=w, actype = actype, activation='relu',pooling=setting.conv_pooling, name='%s_%s%d' % (name, 'EdgeConv', layer_idx))
        
        if mask is not None:
            fts = tf.multiply(fts, mask)

        pool = tf.reduce_mean(fts, axis=1)  # (N, C)
        return pool 
#         if setting.fc_params is not None:
            
#             x = pool
#             for layer_idx, layer_param in enumerate(setting.fc_params):
#                 units, drop_rate = layer_param
#                 x = keras.layers.Dense(units, activation='relu')(x)
#                 if drop_rate is not None and drop_rate > 0:
#                     x = keras.layers.Dropout(drop_rate)(x)
#             out = keras.layers.Dense(setting.num_class, activation='softmax')(x)
# #             out = x
#             return out  # (N, num_classes)
#         else:
#             return pool
        


class _DotDict:
    pass
class _DotDict2:
    pass
class _DotDict3:
    pass
##=================================== MCGCNN =======================================================
class MCECBlock(): #Muti-EdgeConv Block tf.keras.Model
    def __init__(self, points, ptw, features, num_points, K, channels,  w=-1, actype = 'max',with_bn=True, activation='relu', pooling='average', name='edgeconv'):
#         super(MCECBlock, self).__init__(name='')
        ch1, ch2, ch3 = channels, channels, channels
        pts1, pts2, pts3 = tf.identity(points), tf.identity(points), tf.identity(points)
        fts1, fts2, fts3 = tf.identity(features), tf.identity(features), tf.identity(features)
        k1, k2, k3 = K//3, K//3*2, K
        nump1, nump2, nump3 = num_points, num_points, num_points 
        self.econv2a = edge_conv(pts1, ptw, fts1, nump1, k1, ch1, w= w, actype = actype ,with_bn=with_bn, activation=activation, pooling= pooling, name= 'p1_'+name)
        
        self.econv2b = edge_conv(pts2, ptw, fts2, nump2, k2, ch2, w= w, actype = actype,with_bn=with_bn, activation=activation, pooling= pooling, name= 'p2_'+name)
       
        self.econv2c = edge_conv(pts3, ptw, fts3, nump3, k3, ch3, w= w, actype = actype,with_bn=with_bn, activation=activation, pooling= pooling, name= 'p3_'+name)
      
        
    def call(self):
        x1 = self.econv2a
        
        x2 = self.econv2b
        
        x3 = self.econv2c
        
        x = tf.concat([x1,x2,x3], axis=-1)
        x = tf.math.reduce_max(x, axis = -1)
        return x

##==========================================================================================
def get_MCGCNN(num_classes, input_shapes,K=16 , w=-1, actype = 'max' ):
    r"""ParticleNet model from `"ParticleNet: Jet Tagging via Particle Clouds"
    <https://arxiv.org/abs/1902.08570>`_ paper.
    Parameters
    ----------
    num_classes : int
        Number of output classes.
    input_shapes : dict
        The shapes of each input (`points`, `features`, `mask`).
    """
#     MCECBlock = MCECBlock()
    setting = _DotDict()
    setting.num_class = num_classes
    # conv_params: list of tuple in the format (K, (C1, C2, C3))
    setting.conv_params = [
        (K, (64, 64, 64)),
        (K, (128, 128, 128)),
        (K, (256, 256, 256))

        ]
    # conv_pooling: 'average' or 'max'
    setting.conv_pooling = 'average'
    # fc_params: list of tuples in the format (C, drop_rate)
    setting.fc_params = [(256, 0.1)]
    setting.num_points = input_shapes['points'][0]
    
    points = keras.Input(name='points', shape=input_shapes['points'])
    features = keras.Input(name='features', shape=input_shapes['features']) if 'features' in input_shapes else None
    mask = keras.Input(name='mask', shape=input_shapes['mask']) if 'mask' in input_shapes else None
    
    outputs = _particle_net_base(points, features, mask, setting,  w=w, actype = actype, name='MCGCNN1')
    m1 = keras.Model(inputs=[points, features, mask], outputs=outputs, name='MCGCNN1')
# ---------------------------------------------------------
    setting2 = _DotDict2()
    setting2.num_class = num_classes
    # conv_params: list of tuple in the format (K, (C1, C2, C3))
    setting2.conv_params = [
        (K, (64, 64, 64)),
        (K, (128, 128, 128)),
          (K, (256, 256, 256))
        ]
    # conv_pooling: 'average' or 'max'
    setting2.conv_pooling = 'average'
    # fc_params: list of tuples in the format (C, drop_rate)
    setting2.fc_params = [(256, 0.1)]
    setting2.num_points = input_shapes['points'][0]
    out2 = _particle_net_base(points, features, mask, setting2,  w=w+1, actype = actype, name='MCGCNN2')
    m2 = keras.Model(inputs=[points, features, mask], outputs=out2, name='MCGCNN2')
    
    setting3 = _DotDict3()
    setting3.num_class = num_classes
    # conv_params: list of tuple in the format (K, (C1, C2, C3))
    setting3.conv_params = [
        (K, (64, 64, 64)),
        (K, (128, 128, 128)),
          (K, (256, 256, 256))
        ]
    # conv_pooling: 'average' or 'max'
    setting3.conv_pooling = 'average'
    # fc_params: list of tuples in the format (C, drop_rate)
    setting3.fc_params = [(256, 0.1)]
    setting3.num_points = input_shapes['points'][0]
    out3 = _particle_net_base(points, features, mask, setting3,  w=w+1, actype = actype, name='MCGCNN3')
    m3 = keras.Model(inputs=[points, features, mask], outputs=out3, name='MCGCNN3')
    x1 = tf.expand_dims(m1.output, -1)
    x2 = tf.expand_dims(m2.output, -1)
    x3 = tf.expand_dims(m3.output, -1)
#     combined = tf.concat([m1.output, m2.output, m3.output], axis=-1)
    combined = tf.concat([x1, x2, x3], axis=-1)
    combined = tf.expand_dims(combined, -1)
    
    combined = keras.layers.Conv2D(256, kernel_size=(1,1), strides=1, data_format='channels_last',use_bias=True,
                            kernel_initializer='glorot_normal', name='MCGCNN_final_weight')(combined)
    
    combined = keras.layers.BatchNormalization(name='MCGCNN_final_bias')(combined)
            
    combined = keras.layers.Activation('relu', name='MCGCNN_final_relu')(combined)
#     combined = tf.concat([outputs, out2, out3], axis= -1)
    combined =  tf.math.reduce_max(combined, axis = 2)
    combined = keras.layers.Flatten()(combined)
#     combined = tf.squeeze( combined, axis=1)
#     z = tf.math.reduce_max(combined, axis = 1)
    if setting.fc_params is not None:
            
            x = combined 
            for layer_idx, layer_param in enumerate(setting.fc_params):
                units, drop_rate = layer_param
                x = keras.layers.Dense(units, activation='relu')(x)
                if drop_rate is not None and drop_rate > 0:
                    x = keras.layers.Dropout(drop_rate)(x)
    x = keras.layers.Dense(num_classes, activation='relu')(x)
#     z = keras.layers.Dense(128, activation="relu")(combined)
#     z = keras.layers.Dense(num_classes, activation="linear")(x)
    
    
    
    return keras.Model(inputs=[points, features, mask], outputs=x, name='MCGCNN')
##-------------------------------------------------------------------------------------------
def get_particle_net(num_classes, input_shapes,  w=-1, actype = 'max'):
    r"""ParticleNet model from `"ParticleNet: Jet Tagging via Particle Clouds"
    <https://arxiv.org/abs/1902.08570>`_ paper.
    Parameters
    ----------
    num_classes : int
        Number of output classes.
    input_shapes : dict
        The shapes of each input (`points`, `features`, `mask`).
    """
    setting = _DotDict()
    setting.num_class = num_classes
    # conv_params: list of tuple in the format (K, (C1, C2, C3))
    setting.conv_params = [
        (16, (64, 64, 64)),
        (16, (128, 128, 128)),
        (16, (256, 256, 256)),
        ]
    # conv_pooling: 'average' or 'max'
    setting.conv_pooling = 'average'
    # fc_params: list of tuples in the format (C, drop_rate)
    setting.fc_params = [(256, 0.1)]
    setting.num_points = input_shapes['points'][0]

    points = keras.Input(name='points', shape=input_shapes['points'])
    features = keras.Input(name='features', shape=input_shapes['features']) if 'features' in input_shapes else None
    mask = keras.Input(name='mask', shape=input_shapes['mask']) if 'mask' in input_shapes else None
    outputs = _particle_net_base(points, features, mask, setting,  w=-1, actype = actype, name='ParticleNet')

    return keras.Model(inputs=[points], outputs=outputs, name='ParticleNet')


def get_particle_net_lite(num_classes, input_shapes):
    r"""ParticleNet-Lite model from `"ParticleNet: Jet Tagging via Particle Clouds"
    <https://arxiv.org/abs/1902.08570>`_ paper.
    Parameters
    ----------
    num_classes : int
        Number of output classes.
    input_shapes : dict
        The shapes of each input (`points`, `features`, `mask`).
    """
    setting = _DotDict()
    setting.num_class = num_classes
    # conv_params: list of tuple in the format (K, (C1, C2, C3))
    setting.conv_params = [
        (7, (32, 32, 32)),
        (7, (64, 64, 64)),
        ]
    # conv_pooling: 'average' or 'max'
    setting.conv_pooling = 'average'
    # fc_params: list of tuples in the format (C, drop_rate)
    setting.fc_params = [(128, 0.1)]
    setting.num_points = input_shapes['points'][0]

    points = keras.Input(name='points', shape=input_shapes['points'])
    features = keras.Input(name='features', shape=input_shapes['features']) if 'features' in input_shapes else None
    mask = keras.Input(name='mask', shape=input_shapes['mask']) if 'mask' in input_shapes else None
    outputs = _particle_net_base(points, features, mask, setting, name='ParticleNet')

    return keras.Model(inputs=[points, features, mask], outputs=outputs, name='ParticleNet')
