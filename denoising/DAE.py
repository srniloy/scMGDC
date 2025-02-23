import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, Dropout, BatchNormalization, Add
from tensorflow.keras.models import Model
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import scanpy as sc
import numpy as np

def add_sparse_noise(data, dropout_rate=0.1):
    noisy_data = data.copy()
    mask = np.random.binomial(1, dropout_rate, size=data.shape)  # Create dropout mask
    noisy_data[mask == 1] = 0  # Set selected values to zero
    return noisy_data

def autoencoder():
        
    X = adata.X
    X = X.toarray() if hasattr(X, "toarray") else X  # Convert to dense if sparse

    X = sc.pp.log1p(X, copy=True)
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    X_noisy = add_sparse_noise(X)

    input_dim = X.shape[1]

    input_layer = Input(shape=(input_dim,))
    encoded = Dense(512, activation='relu')(input_layer)
    encoded = BatchNormalization()(encoded)
    encoded = Dropout(0.2)(encoded)
    latent_space = Dense(128, activation='relu')(encoded)

    decoded = Dense(512, activation='relu')(latent_space)
    decoded = BatchNormalization()(decoded)
    decoded = Dropout(0.2)(decoded)
    output_layer = Dense(input_dim, activation='sigmoid')(decoded)

    autoencoder = Model(inputs=input_layer, outputs=output_layer)

    def sparse_mse(y_true, y_pred):
        mask = tf.cast(tf.equal(y_true, 0.0), tf.float32)  # Mask for zeros in original data
        nonzero_loss = tf.reduce_mean(tf.square((1 - mask) * (y_true - y_pred)))  # Penalize changes to nonzeros
        zero_loss = tf.reduce_mean(tf.square(mask * y_pred))  # Penalize predicted values where true value is zero
        return nonzero_loss + 0.1 * zero_loss  # Combine losses with weight for zero preservation

    # Compile the model
    autoencoder.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=1e-3),
                        loss=sparse_mse,
                        metrics=['mae'])

    # Step 4: Train the model
    X_train, X_val = train_test_split(X, test_size=0.2, random_state=42)
    X_noisy_train, X_noisy_val = train_test_split(X_noisy, test_size=0.2, random_state=42)

    history = autoencoder.fit(
        X_noisy_train, X_train,
        validation_data=(X_noisy_val, X_val),
        epochs=200,
        batch_size=256,
        shuffle=True,
        callbacks=[
            tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=15, restore_best_weights=True),
            tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5)
        ]
    )

    # Step 5: Denoise the data
    X_denoised = autoencoder.predict(X_noisy)

    # Save the denoised output
    adata.X = scaler.inverse_transform(X_denoised)  # Inverse scaling
    
    return adata