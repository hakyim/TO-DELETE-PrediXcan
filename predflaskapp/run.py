#!flask/bin/python
from app import app
import ssl

context =ssl.SSLContext(ssl.PROTOCOL_TLSv1_2)
context.load_cert_chain('yourserver.crt', 'yourserver.key')

app.run(debug=True,ssl_context=context)